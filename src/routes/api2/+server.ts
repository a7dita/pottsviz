import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { exec } from 'node:child_process';
import fs from 'node:fs/promises';
import util from 'node:util';

const execPromisified = util.promisify(exec);

// this helper function uses the promosified exec method to run the sh command.
// I do not understand fully.
const runCommand = async (command: string) => {
    const { stdout, stderr } = await execPromisified(command);
    return stdout;
};

// upon receiving a POST request, this api runs the command and send the result to the frontend.
export const POST: RequestHandler = async ({ request }) => {
    const { command } = await request.json();
    const result = await runCommand(command);
    console.log(result)
    return new Response(JSON.stringify(result));

};
// FIXME we need to edit automate.py to get the PID of the c++ process as the 'result',
// so that we can write another POST method to pkill the process upon 'stop' button pressing.
