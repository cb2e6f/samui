#!/usr/bin/env python
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

#%%
# You need to download the data from https://info.vizgen.com/mouse-brain-data
# which, for each slice, contains the following files:
#   cell_boundaries/
#   cell_by_gene_S{n}R{n}.csv
#   cell_metadata_S{n}R{n}.csv
#   detected_transcripts_S1R1.csv
#   images/
#   |- micron_to_mosaic_pixel_transform.csv
#   |- mosaic_DAPI_z{n}.tif
#
# where {n} is a number.
# This script expects these files in a single directory (use one DAPI file).
# Note that we need to convert the uint16 format in of the DAPI image to uint8
# for the compression.
# This requires a huge amount of memory since each image is 10.2 GB.
# To download data from Google Cloud in the command line, use [gsutil](https://cloud.google.com/storage/docs/downloading-objects).

# sample_dir = Path("/Users/chaichontat/Downloads/trs")
sample_dir = sys.argv[1]
# dapi = (
#     sample_dir / "datasets-mouse_brain_map-BrainReceptorShowcase-Slice1-Replicate1-images-mosaic_DAPI_z0.tif"
# )

print(sample_dir)

metadata_file = open(os.path.join(sample_dir, "cell_metadata.csv"), 'r')
by_gene_file = open(os.path.join(sample_dir, "cell_by_gene.csv"), 'r')

entityid_map_file = open(sample_dir +  "/" + "entityid_map.csv", 'w')
mapped_metadata_file = open(sample_dir +  "/" + "mapped_cell_metadata.csv", 'w')
mapped_by_gene_file = open(sample_dir +  "/" + "mapped_cell_by_gene.csv", 'w')

count = 0

for metadata_line in metadata_file:

    by_gene_line = by_gene_file.readline()

    if count == 0:
        mapped_metadata_file.write(metadata_line)
        mapped_by_gene_file.write(by_gene_line)
        count = 1
        continue

    metadata_line_split = metadata_line.split(',')
    by_gene_split = by_gene_line.split(',')

    if metadata_line_split[0] != by_gene_split[0]:
        print("entity id mismatch - this isn't going to work!")

    entityid_map_file.write(str(count) + "," + metadata_line_split[0] + "\n")

    metadata_line_split[0] = str(count)
    mapped_metadata_file.write(",".join(metadata_line_split))

    by_gene_split[0] = str(count)
    mapped_by_gene_file.write(",".join(by_gene_split))

    count = count + 1

mapped_by_gene_file.close()
mapped_metadata_file.close()

by_gene_file.close()
metadata_file.close()
entityid_map_file.close()

#%% Coords
coords = remove_dupes(
    pd.read_csv(sample_dir +  "/" + "mapped_cell_metadata.csv", index_col=0, dtype=np.float32).rename(
        columns={"center_x": "x", "center_y": "y"}
    )[["x", "y"]]
)
print(coords)
coords.index = coords.index.map("{:.0f}".format)

feat = remove_dupes(
    pd.read_csv(
        sample_dir + "/" + "mapped_cell_by_gene.csv",
        index_col=0,
        dtype=np.float32,
    )
).apply(lambda x: np.log2(x + 1), raw=True, axis=0)
print(feat)
feat.index = feat.index.map("{:.0f}".format)

shutil.rmtree(sample_dir + "/" + "out_image", ignore_errors=True)
# This creates a new sample that has only spots.
s = (
    Sample(name="sample", path=sample_dir + "/" + "out")
    # You may need to adjust the size to make the spot size more reasonable.
    .add_coords(coords, name="cellCoords", mPerPx=1e-6, size=1e-5).add_chunked_feature(
        feat, name="cells", coordName="cellCoords", unit="Log counts", sparse=True
    )
    .set_default_feature(group="cells", feature="TraesCS7D02G261600") # Example: this needs to match a feature that exists in your sample.
    .write()
)

# You can drag the result in the out folder to Samui to see your spots.

# %%
# Next, we will add the image to the sample.
# This requires alignment from the spot coordinates to the image coordinates.

# Affine matrix
scale = np.loadtxt(sample_dir + "/" + "images" + "/" + "micron_to_mosaic_pixel_transform.csv")

# Inverse transform
affine = ~Affine(*scale[:2].flatten() * 1e6)
dapi = sample_dir + "/" + "images" + "/" + "mosaic_DAPI_z3.tif"


# Make another sample that has an image.
shutil.rmtree(sample_dir + "/" + "out_image", ignore_errors=True)
shutil.copytree(sample_dir + "/" + "out", sample_dir + "/" + "out_image")
s.set_path(Path(sample_dir + "/" + "out_image"))

#%%
s.add_image(dapi, channels=["DAPI"], scale=affine.a, translate=(affine.c * 1e-6, affine.f * 1e-6)).write()

# %%
