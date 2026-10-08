import { cpSync, readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

// Native programs are not JavaScript imports: explicitly package them in each
// generated function instead of relying on automatic dependency tracing.
function bundle(directory) {
  for (const entry of readdirSync(directory, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue;
    const path = join(directory, entry.name);
    if (entry.name.endsWith('.func')) {
      cpSync('.simulators', join(path, '.simulators'), { recursive: true });
      const configPath = join(path, '.vc-config.json');
      const config = JSON.parse(readFileSync(configPath, 'utf8'));
      config.supportsCancellation = true;
      writeFileSync(configPath, JSON.stringify(config, null, 2));
    } else bundle(path);
  }
}
bundle('.vercel/output/functions');
