from pathlib import Path
from typing import Callable, Literal, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import validator
from scipy.sparse import csc_matrix, csr_matrix
from typing_extensions import Self

from .utils import ReadonlyModel, Url, Writable, concat

FeatureType = Literal["categorical", "quantitative", "singular"]


class Coord(ReadonlyModel):
    x: int
    y: int


class CoordId(Coord):
    id: str


class CoordParams(Writable):
    """
    name: Name of the overlay
    type: Type of the overlay 'single', 'multi'
    shape: Shape of the overlay (currently only circle)
    url: points to a json file with the coordinates of the overlay
        in {x: number, y: number, id?: string}[].
        or
        ChunkedHeader
    mPerPx: micrometers per pixel
    size: size of the overlay in micrometers
    """

    name: str
    shape: Literal["circle"]
    url: Url
    mPerPx: float | None = None
    size: float | None = None


class PlainCSVParams(Writable):
    type: Literal["plainCSV"] = "plainCSV"
    name: str
    url: Url
    dataType: FeatureType = "quantitative"
    coordName: str | None = None
    unit: str | None = None
    size: float | None = None


class ChunkedCSVParams(ReadonlyModel):
    type: Literal["chunkedCSV"] = "chunkedCSV"
    name: str
    url: Url
    headerUrl: Url | None = None
    dataType: FeatureType = "quantitative"
    unit: str | None = None

    @validator("headerUrl", always=True, pre=False)
    def check_headerUrl(cls, v: Url | None, values: dict[str, str | None]):
        if v is None:
            path = Path(cast(Url, values["url"]).url)
            return Url(url=path.with_suffix(".json").as_posix())
        return v

    def write(self, f: Callable[[Path], None], header: Callable[[Path], None] | None = None) -> Self:
        f(Path(self.url.url))
        if header and self.headerUrl:
            header(Path(self.headerUrl.url))
        return self


class ChunkedCSVHeader(ReadonlyModel):
    names: list[str] | None = None
    ptr: list[int]
    length: int
    activeDefault: str | None = None
    sparseMode: Literal["record", "array"] | None = None
    weights: list[float] | None = None
    coordName: str | None = None


class FeatureAndGroup(ReadonlyModel):
    feature: str
    group: str | None = None


FeatureParams = ChunkedCSVParams | PlainCSVParams


def get_compressed_genes(
    arr: npt.ArrayLike | csc_matrix | csr_matrix, names: list[str], *, coordName: str, mode: Literal["csr", "csc"] = "csc"  # type: ignore
) -> tuple[ChunkedCSVHeader, bytearray]:
    if mode == "csr":
        cs = csr_matrix(arr)  # csR
    elif mode == "csc":
        cs = csc_matrix(arr)  # csC
    else:
        raise ValueError("Invalid mode")

    indices = cs.indices.astype(int)
    indptr = cs.indptr.astype(int)
    data = cs.data

    objs = []

    for i in range(len(indptr) - 1):
        if (indices[indptr[i] : indptr[i + 1]]).size == 0:
            objs.append(None)
        else:
            objs.append(
                pd.DataFrame(
                    {
                        "index": indices[indptr[i] : indptr[i + 1]].tolist(),
                        "value": [round(x, 3) for x in data[indptr[i] : indptr[i + 1]].tolist()],
                    }
                )
            )

    ptr, outbytes = concat(objs, lambda x: x.to_csv(index=False).encode())
    match mode:
        case "csr":
            length = cs.shape[1]
        case "csc":
            length = cs.shape[0]
        case _:
            raise ValueError("Unknown mode")

    return (
        ChunkedCSVHeader(
            names=names,
            weights=list(
                np.array(cs.sum(axis=0)).flatten().round(2)
                if mode == "csc"
                else np.array(cs.sum(axis=1)).flatten().round(2)
            ),
            ptr=ptr.tolist(),
            length=length,
            sparseMode="array" if mode == "csc" else "record",
            coordName=coordName,
        ),
        outbytes,
    )
