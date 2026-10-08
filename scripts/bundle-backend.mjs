import { cpSync, readdirSync } from 'node:fs';
import { join } from 'node:path';

// Native programs are not JavaScript imports: explicitly package them in each
// generated function instead of relying on automatic dependency tracing.
function bundle(directory) {
  for (const entry of readdirSync(directory, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue;
    const path = join(directory, entry.name);
    if (entry.name.endsWith('.func')) {
      cpSync('.simulators', join(path, '.simulators'), { recursive: true });
    } else bundle(path);
  }
}
bundle('.vercel/output/functions');
