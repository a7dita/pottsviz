<script lang="ts">
  import type { topologies } from '$lib/topologies';
  import katex from 'katex';
  export let topology: typeof topologies[number];
  const math = (value: string) => katex.renderToString(value, { throwOnError: false });
  const definitions = [
    { key: 'redundancy_log1p', label: 'Redundancy', symbol: 'R_4',
      formula: String.raw`R_4=\log(1+N_4),\quad N_4=\sum_{f<g}\binom{(BB^\top)_{fg}}{2}`,
      text: 'Counts four-edge cycles: two frontal memories sharing two posterior neighbours. The score is log(1 + cycle count). Higher values mean more shared associations.' },
    { key: 'q_star', label: 'Modularity', symbol: 'Q^*',
      formula: String.raw`Q^*\approx\max_c\frac{1}{49k}\sum_{f,p}\left(B_{fp}-\frac{k}{49}\right)\delta_{c_f,c_p}`,
      text: 'Barber bipartite modularity, estimated by the pilot’s fixed-restart greedy search. Higher values mean stronger communities compared with a degree-preserving null. T00 is outside this degree-seven comparison; its manifest leaves Q* undefined.' },
    { key: 'mixing_gap', label: 'Mixing gap', symbol: '\\gamma',
      formula: String.raw`\gamma=1-\sigma_2(B/k)^2`,
      text: 'Spectral gap of a two-step memory walk. σ₂ is the second singular value of the normalized association matrix. A larger gap means faster mixing; zero indicates disconnected components.' }
  ];
</script>

<div class="topology-view">
  <p class="caption">Frontal memories ↓ · posterior memories →</p>
  <svg viewBox="0 0 310 310" role="img" aria-label="{topology.id} association matrix: 49 frontal by 49 posterior memories">
    <rect x="12" y="12" width="294" height="294" fill="#f3f0fa" />
    {#each topology.matrix as row, f}
      {#each row as edge, p}
        {#if edge}<rect x={12+p*6} y={12+f*6} width="5.5" height="5.5" fill="#6d28d9"><title>Frontal {f} ↔ posterior {p} · weight {topology.edge_weight_after_loading.toFixed(4)}</title></rect>{/if}
      {/each}
    {/each}
    {#each [0, 7, 14, 21, 28, 35, 42] as index}
      <text x={12+index*6} y="9" font-size="8" fill="#6b7280">{index}</text>
      <text x="10" y={17+index*6} text-anchor="end" font-size="8" fill="#6b7280">{index}</text>
    {/each}
  </svg>
  <p class="caption">49 + 49 memories · {topology.degree} link{topology.degree === 1 ? '' : 's'} per memory · weight 1/{topology.degree}</p>
  <div class="scores">
    {#each definitions as definition, i}
      <div class="score">
        <button type="button" aria-describedby="score-help-{i}">
          <span>{definition.label} {@html math(definition.symbol)} ⓘ</span>
          <strong>{topology.metrics[definition.key as keyof typeof topology.metrics]?.toFixed(3) ?? '—'}</strong>
        </button>
        <div class="help" id="score-help-{i}" role="tooltip">
          <p>{definition.text}</p><div class="formula">{@html math(definition.formula)}</div>
          {#if i === 0}<p>Four-edge cycles N₄ = {topology.metrics.n4}.</p>{/if}
          <p>B is the binary frontal × posterior matrix; k is its degree.</p>
        </div>
      </div>
    {/each}
  </div>
  <p class="caption">Hover or focus a score for its definition and equation.</p>
</div>

<style>
  svg { width: 100%; }
  .caption { font-size: .75rem; color: #6b7280; margin: .5rem 0; }
  .scores { display: grid; gap: .5rem; }
  .score { position: relative; }
  button { display: flex; justify-content: space-between; gap: .5rem; width: 100%; text-align: left; background: #f3f0fa; padding: .6rem; border-radius: .4rem; font-size: .85rem; }
  button:focus-visible { outline: 2px solid #7c3aed; }
  .help { display: none; position: absolute; z-index: 5; left: 0; top: 100%; width: min(460px, 80vw); background: white; padding: 1rem; border: 1px solid #ddd6fe; box-shadow: 0 6px 20px #0002; border-radius: .5rem; font-size: .8rem; }
  .score:hover .help, .score:focus-within .help { display: block; }
  .help p { margin: .5rem 0; }
  .formula { overflow-x: auto; padding: .5rem 0; }
</style>
