import React, { useState, useEffect } from 'react';
import { Atom, Zap, Database, Users, ArrowRight, Play, Check } from 'lucide-react';
import SolutionsSection from '../components/SolutionsSection';

export default function QuantumChemistryHomepage() {
  const [scrollY, setScrollY] = useState(0);
  const [activeFeature, setActiveFeature] = useState(0);

  useEffect(() => {
    // Remove default margins from body
    document.body.style.margin = '0';
    document.body.style.padding = '0';
    document.documentElement.style.margin = '0';
    document.documentElement.style.padding = '0';
    
    const handleScroll = () => setScrollY(window.scrollY);
    window.addEventListener('scroll', handleScroll);
    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  useEffect(() => {
    const interval = setInterval(() => {
      setActiveFeature((prev) => (prev + 1) % 4);
    }, 3000);
    return () => clearInterval(interval);
  }, []);

  const features = [
    {
      icon: <Atom className="w-8 h-8" />,
      title: "Molecular Dynamics",
      description: "Simulate complex molecular interactions with precision and speed"
    },
    {
      icon: <Zap className="w-8 h-8" />,
      title: "Fast Computing",
      description: "GPU-accelerated calculations for real-time results"
    },
    {
      icon: <Database className="w-8 h-8" />,
      title: "Data Management",
      description: "Organize and analyze simulation data effortlessly"
    },
    {
      icon: <Users className="w-8 h-8" />,
      title: "Collaboration",
      description: "Share simulations and insights with your research team"
    }
  ];

  const useCases = [
    "Drug Discovery & Development",
    "Materials Science Research",
    "Catalysis Optimization",
    "Protein Structure Analysis"
  ];

  return (
    <div className="min-h-screen min-w-screen bg-gradient-to-br from-slate-900 via-purple-900 to-slate-900 text-white overflow-hidden">
      {/* Animated Background */}
      <div className="fixed inset-0 overflow-hidden pointer-events-none">
        <div 
          className="absolute w-96 h-96 bg-purple-500/20 rounded-full blur-3xl"
          style={{
            top: '10%',
            left: '20%',
            transform: `translate(${scrollY * 0.1}px, ${scrollY * 0.15}px)`
          }}
        />
        <div 
          className="absolute w-96 h-96 bg-blue-500/20 rounded-full blur-3xl"
          style={{
            bottom: '10%',
            right: '20%',
            transform: `translate(-${scrollY * 0.08}px, -${scrollY * 0.12}px)`
          }}
        />
      </div>

      {/* Navigation */}
      <nav className="relative z-10 px-6 py-6 flex items-center justify-between max-w-7xl mx-auto">
        <div className="flex items-center space-x-2">
          <Atom className="w-8 h-8 text-purple-400" />
          <span className="text-2xl font-bold bg-gradient-to-r from-purple-400 to-blue-400 bg-clip-text text-transparent">
            QuantumMLD
          </span>
        </div>
        <div className="hidden md:flex space-x-8 text-sm">
          <div className="hidden md:flex space-x-8 text-sm">
            <a href="/" className="hover:text-purple-400 transition">Home</a>
            <a href="/plots" className="hover:text-purple-400 transition">Plot Generator</a>
            <a href="#features" className="hover:text-purple-400 transition">Features</a>
            <a href="#docs" className="hover:text-purple-400 transition">Docs</a>
          </div>
        </div>
        <div className="flex space-x-4">
          <button className="px-4 py-2 text-sm hover:text-purple-400 transition">
            Sign In
          </button>
          <button className="px-6 py-2 bg-gradient-to-r from-purple-600 to-blue-600 rounded-lg text-sm font-medium hover:shadow-lg hover:shadow-purple-500/50 transition">
            Get Started
          </button>
        </div>
      </nav>

      {/* Hero Section */}
      <section className="relative z-10 px-6 pt-20 pb-32 max-w-7xl mx-auto">
        <div className="text-center space-y-8">
          <div className="inline-block">
            <div className="px-4 py-2 bg-purple-500/20 rounded-full border border-purple-500/30 text-sm">
              ✨ Go powered developemet
            </div>
          </div>
          
          <h1 className="text-6xl md:text-7xl font-bold leading-tight">
            Quantum Chemistry
            <br />
            <span className="bg-gradient-to-r from-purple-400 via-pink-400 to-blue-400 bg-clip-text text-transparent">
              Reimagined
            </span>
          </h1>
          
          <p className="text-xl text-gray-300 max-w-2xl mx-auto">
            Run sophisticated quantum simulations in the cloud. Accelerate your research 
            with our powerful, intuitive platform built for modern chemistry.
          </p>

          <div className="flex flex-col sm:flex-row items-center justify-center gap-4 pt-8">
            <button className="group px-8 py-4 bg-gradient-to-r from-purple-600 to-blue-600 rounded-lg font-medium text-lg hover:shadow-2xl hover:shadow-purple-500/50 transition flex items-center space-x-2">
              <Play className="w-5 h-5" />
              <span>Try Yourself</span>
              <ArrowRight className="w-5 h-5 group-hover:translate-x-1 transition" />
            </button>
            <button className="px-8 py-4 border border-purple-500/50 rounded-lg font-medium text-lg hover:bg-purple-500/10 transition">
              Watch Demo
            </button>
          </div>

          {/* Stats */}
          <div className="grid grid-cols-3 gap-8 max-w-3xl mx-auto pt-16">
            <div>
              <div className="text-4xl font-bold text-purple-400">1google</div>
              <div className="text-sm text-gray-400 mt-1">Simulations Run</div>
            </div>
            <div>
              <div className="text-4xl font-bold text-blue-400">4</div>
              <div className="text-sm text-gray-400 mt-1">Developers</div>
            </div>
            <div>
              <div className="text-4xl font-bold text-pink-400">100%</div>
              <div className="text-sm text-gray-400 mt-1">1D Simulations</div>
            </div>
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section id="features" className="relative z-10 px-6 py-20 bg-black/20 backdrop-blur-sm">
        <div className="max-w-7xl mx-auto">
          <div className="text-center mb-16">
            <h2 className="text-4xl font-bold mb-4">Powerful Features</h2>
            <p className="text-gray-400 text-lg">Everything you need for cutting-edge quantum simulations</p>
          </div>

          <div className="grid md:grid-cols-2 lg:grid-cols-4 gap-6">
            {features.map((feature, index) => (
              <div
                key={index}
                className={`p-6 rounded-2xl border transition-all duration-500 cursor-pointer ${
                  activeFeature === index
                    ? 'bg-gradient-to-br from-purple-600/20 to-blue-600/20 border-purple-500/50 scale-105'
                    : 'bg-white/5 border-white/10 hover:border-purple-500/30'
                }`}
                onMouseEnter={() => setActiveFeature(index)}
              >
                <div className="text-purple-400 mb-4">{feature.icon}</div>
                <h3 className="text-xl font-semibold mb-2">{feature.title}</h3>
                <p className="text-gray-400 text-sm">{feature.description}</p>
              </div>
            ))}
          </div>
        </div>
      </section>

      {/* Use Cases */}
      <SolutionsSection />

      {/* Footer */}
      <footer className="relative z-10 px-6 py-12 border-t border-white/10">
        <div className="max-w-7xl mx-auto text-center text-gray-400 text-sm">
          <p>© 2025 QuantumMLD. Brewed to Compute.</p>
        </div>
      </footer>
    </div>
  );
}