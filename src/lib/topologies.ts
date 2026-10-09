import manifest from './topologies/manifest.json';

const csvFiles = import.meta.glob('./topologies/*.csv', { eager: true, query: '?raw', import: 'default' });
export const topologies = manifest.topologies.map(item => ({
  ...item,
  matrix: String(csvFiles[`./topologies/${item.file}`]).trim().split('\n').map(row => row.split(',').map(Number))
}));

export type TopologyFamily = 'original' | 'random' | 'modular' | 'shared-target';
export function topologyForFamily(family: TopologyFamily, parameter: number = 0) {
  const id = family === 'modular' ? (parameter === 7 ? 'ms7' : `m${parameter}`)
    : family === 'shared-target' ? (parameter === 7 ? 'ms7' : `s${parameter}`) : family;
  const topology = topologies.find(item => item.id === id);
  if (!topology) throw new Error(`Unknown topology: ${id}`);
  return topology;
}
