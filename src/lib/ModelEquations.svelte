<script lang="ts">
  import katex from 'katex';
  import 'katex/dist/katex.min.css';
  import { modelEquations } from './model-equations';

  export let model: keyof typeof modelEquations;
  $: equations = modelEquations[model];
</script>

<section class="equations" aria-label="Key model equations">
  <h2>Model equations</h2>
  <div class="rows">
    {#each equations.rows as row}
      <div class="row">
        <p class="label">{row.label}</p>
        <div class="math">{@html katex.renderToString(row.latex, { displayMode: true, output: 'mathml', throwOnError: true })}</div>
      </div>
    {/each}
  </div>
  <p class="note">{equations.note} The simulator applies unit-wise Euler updates with one time step per sweep.</p>
</section>

<style>
  .equations { min-width: 0; width: 100%; border: 1px solid #e9d5ff; border-radius: .6rem; background: #fdfbff; padding: 1rem 1.25rem; box-sizing: border-box; text-align: left; }
  h2 { font-size: 1rem; font-weight: 600; color: #6d28d9; margin-bottom: .7rem; }
  .rows { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 390px), 1fr)); gap: .6rem 1.4rem; }
  .row { min-width: 0; }
  .label { font-size: .72rem; color: #7c718e; margin: 0 0 .2rem; }
  .math { overflow-x: auto; overflow-y: hidden; padding: .25rem 0 .65rem; color: #4b3568; }
  .math :global(math) { font-family: 'STIX Two Math', 'Cambria Math', math; font-size: 1rem; white-space: nowrap; }
  .note { font-size: .75rem; line-height: 1.5; color: #756d82; margin-top: .45rem; }
  @media (max-width: 700px) { .equations { padding: .8rem; } .math :global(math) { font-size: .9rem; } }
</style>
