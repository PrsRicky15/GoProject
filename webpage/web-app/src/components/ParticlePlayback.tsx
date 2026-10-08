import { useEffect, useState } from 'react';
import type { ParticleFrame } from '../lib/simulation';

export default function ParticlePlayback({ frames, box }: { frames: ParticleFrame[]; box: number }) {
  const [index, setIndex] = useState(0);
  const [playing, setPlaying] = useState(false);
  useEffect(() => {
    if (!playing) return;
    const timer = window.setInterval(() => setIndex(i => (i + 1) % frames.length), 80);
    return () => window.clearInterval(timer);
  }, [playing, frames.length]);
  const frame = frames[index];
  if (!frame) return null;
  return <section className="panel particle-panel"><div className="section-heading"><h3>Particle playback</h3><span className="tag">Time {frame.time.toFixed(3)}</span></div>
    <svg viewBox={`0 0 ${box} ${box}`} className="particle-box" role="img" aria-label={`${frame.x.length} particles at time ${frame.time.toFixed(3)} in a periodic box`}><rect width={box} height={box} fill="#f6f3fc" />{frame.x.map((x, i) => <circle key={i} cx={x} cy={box - frame.y[i]} r={.18} fill={i === 0 ? '#f59b46' : '#7555d8'} />)}</svg>
    <div className="playback-controls"><button type="button" onClick={() => setPlaying(!playing)}>{playing ? 'Pause' : 'Play'}</button><input aria-label="Simulation frame" type="range" min="0" max={frames.length - 1} value={index} onChange={e => { setPlaying(false); setIndex(Number(e.target.value)); }} /><span>{index + 1}/{frames.length}</span></div>
    <p className="muted">Opposite edges are connected. Playback uses the sampled Go trajectory; it does not integrate the equations in the browser.</p>
  </section>;
}
