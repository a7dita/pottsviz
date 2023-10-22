import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { exec } from 'node:child_process';
import fs from 'node:fs/promises';
import util from 'node:util';

const execPromisified = util.promisify(exec);

const runCommand = async (command: string) => {
    const { stdout, stderr } = await execPromisified(command);
    return stdout;
};

export const POST: RequestHandler = async ({ request }) => {
    const { command } = await request.json();
    const result = await runCommand(command);
    console.log(result)
    // return new Response(result);
    return new Response(JSON.stringify(result));

};
