import manifest from './topologies/manifest.json';

const csvFiles = import.meta.glob('./topologies/*.csv', { eager: true, query: '?raw', import: 'default' });
export const topologies = manifest.topologies.map(item => ({
  ...item,
  matrix: String(csvFiles[`./topologies/${item.file}`]).trim().split('\n').map(row => row.split(',').map(Number))
}));
