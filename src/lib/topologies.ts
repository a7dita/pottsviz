import manifest from './topologies/manifest.json';
import displayOrders from './topology-display-order.json';

const csvFiles = import.meta.glob('./topologies/*.csv', { eager: true, query: '?raw', import: 'default' });
export const topologies = manifest.topologies.map(item => ({
  ...item,
  displayOrder: displayOrders[item.id as keyof typeof displayOrders],
  matrix: String(csvFiles[`./topologies/${item.file}`]).trim().split('\n').map(row => row.split(',').map(Number))
}));
