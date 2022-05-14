import sampleList from '$lib/data/meh';
import { get, writable, type Writable } from 'svelte/store';
import type { NameWithFeature } from './data/features';
import type { Sample } from './data/sample';

export type Idx = { idx: number; source: string };

export type State = {
  lockedIdx: Idx;
  currIdx: Idx;
  readonly locked: boolean;
};

export const store: Writable<State> = writable({
  lockedIdx: { idx: -1, source: 'scatter' },
  currIdx: { idx: 0, source: 'scatter' },
  get locked() {
    return this.lockedIdx.idx !== -1;
  }
});

export type HoverName<T> = {
  hover: T | null;
  selected: T | null;
  get active(): T | null;
};

export type NameWithFeatures = {
  feature?: string;
  names: string[];
};

export function genHoverName<T>({ hover, selected }: { hover?: T; selected?: T }): HoverName<T> {
  return {
    hover: hover ?? null,
    selected: selected ?? null,
    get active() {
      if (this.hover) return this.hover;
      return this.selected;
    }
  };
}

export const samples: Writable<Record<string, Sample>> = writable({});

type OverlayName = string;
export const activeFeatures: Writable<Record<OverlayName, NameWithFeature>> = writable({});
export const activeOverlay: Writable<string> = writable('');

// Changed in mapTile and byod.
// Automatically hydrates on change.
export const activeSample: Writable<string> = writable('');
activeSample.subscribe((n) => {
  const s = get(samples)[n];
  if (!s) return;
  if (!s.hydrated) s.hydrate().catch(console.error);
});

export const activeMap: Writable<number> = writable(0);
export const mapList: Writable<number[]> = writable([]);

export function preload() {
  if (sampleList) {
    sampleList
      .then((ss) => {
        const obj = Object.fromEntries(ss.map((s) => [s.name, s]));
        samples.set({ ...obj, ...get(samples) });
      })
      .catch(console.error);
  }
}
