<script lang="ts">
  import type { Sample } from './simulation-client';
  export let samples: Sample[] = [];
  export let region: 'posterior' | 'frontal' = 'posterior';
  const width = 700, height = 440;
  const left = 52, right = 18, top = 18, bottom = 48;
  $: maxTime = Math.max(1, samples[samples.length - 1]?.time ?? 0);
  $: memories = samples[0]?.[region]?.length ?? 0;
  $: paths = Array.from({ length: memories }, (_, memory) => samples.map((sample, i) => {
    const value = sample[region]?.[memory];
    if (value === undefined || !Number.isFinite(value)) return '';
    const x = left + sample.time / maxTime * (width - left - right);
    const y = top + (1 - value) / 1.1 * (height - top - bottom);
    return `${i ? 'L' : 'M'}${x.toFixed(2)},${y.toFixed(2)}`;
  }).join(' '));
</script>

<svg viewBox="0 0 {width} {height}" role="img" aria-label="Live {region} memory overlaps" style="width:100%;max-height:440px">
  <defs><clipPath id="clip-{region}"><rect x={left} y={top} width={width-left-right} height={height-top-bottom}/></clipPath></defs>
  {#each [0, 0.25, 0.5, 0.75, 1] as value}
    <line x1={left} x2={width-right} y1={top+(1-value)/1.1*(height-top-bottom)} y2={top+(1-value)/1.1*(height-top-bottom)} stroke="#e5e7eb"/>
    <text x={left-8} y={top+(1-value)/1.1*(height-top-bottom)+4} text-anchor="end" font-size="12" fill="#4b5563">{value}</text>
  {/each}
  <g clip-path="url(#clip-{region})">
    {#each paths as path, i}<path d={path} stroke="hsl({i*137.508%360},70%,42%)" stroke-width="1.3" fill="none" />{/each}
  </g>
  <line x1={left} x2={width-right} y1={height-bottom} y2={height-bottom} stroke="#6b7280"/>
  <text x={left} y={height-bottom+20} font-size="12">0</text>
  <text x={width-right} y={height-bottom+20} text-anchor="end" font-size="12">{maxTime}</text>
  <text x={width/2} y={height-6} text-anchor="middle" font-size="13">Time (network sweeps)</text>
  <text transform="translate(14,{height/2}) rotate(-90)" text-anchor="middle" font-size="13">Memory overlap</text>
  {#if !samples.length}<text x={width/2} y={height/2} text-anchor="middle" fill="#6b7280">Start a simulation to see the memory overlaps.</text>{/if}
</svg>
