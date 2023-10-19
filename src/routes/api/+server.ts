import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { exec } from 'node:child_process';
import fs from 'node:fs/promises';
import util from 'node:util';

// const execPromisified = util.promisify(exec);

// const runCommand = async (command: string) => {
// 	const { stdout, stderr } = await execPromisified(command);
// 	return stdout;
// };

const readOutput = async () => {
	try {
		const output = await fs.readFile('data/mall_cue0_seed1_d', 'utf8');
		return output;
	} catch (err) {
		return err
	}
}

export const POST: RequestHandler = async ({ request }) => {
	// const { command } = await request.json();
	// const result = await runCommand(command);
	// console.log(result)
	const output = await readOutput();
	return new Response(output);

};
