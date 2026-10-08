import { lazy, Suspense, useEffect, useState } from 'react'
import Homepage from './pages/HomePage'
import './index.css'

const Workbench = lazy(() => import('./pages/Workbench'))
const Calculations = lazy(() => import('./pages/Calculations'))

function App() {
  const [hash, setHash] = useState(() => window.location.hash)
  useEffect(() => {
    const update = () => setHash(window.location.hash)
    window.addEventListener('hashchange', update)
    return () => window.removeEventListener('hashchange', update)
  }, [])
  return hash === '#tools' || hash === '#calculations'
    ? <Suspense fallback={<p role="status" style={{ padding: 32 }}>Opening scientific workspace…</p>}>{hash === '#calculations' ? <Calculations /> : <Workbench />}</Suspense>
    : <Homepage />
}

export default App
