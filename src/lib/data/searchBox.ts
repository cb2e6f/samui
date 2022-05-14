import type { NameWithFeatures } from '../store';
import { ChunkedJSON, PlainJSON, type Data } from './features';

export function updateNames(
  features: Record<string, Data>,
  filterOverlay: string
): NameWithFeatures[] {
  if (!features) return [];
  const out: NameWithFeatures[] = [{ feature: undefined, names: [] }];
  for (const [name, f] of Object.entries(features)) {
    if (f.overlay !== filterOverlay) continue;

    if (f instanceof PlainJSON) {
      out[0].names.push(name);
    } else if (f instanceof ChunkedJSON) {
      if (f.header!.names) {
        out.push({ feature: name, names: Object.keys(f.header!.names) });
      }
    } else {
      throw new Error('Unknown feature type');
    }
  }
  return out;
}
