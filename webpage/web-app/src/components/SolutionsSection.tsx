import React from 'react';
import { Check } from 'lucide-react';

export default function SolutionsSection() {
  const useCases = [
    "Drug Discovery & Development",
    "Materials Science Research",
    "Catalysis Optimization",
    "Protein Structure Analysis"
  ];

  return (
    <section id="solutions" className="relative z-10 px-6 py-20 max-w-7xl mx-auto">
      <div className="text-center mb-16">
        <h2 className="text-4xl font-bold mb-4">Trusted Across Industries</h2>
        <p className="text-gray-400 text-lg">From pharmaceuticals to advanced materials</p>
      </div>

      <div className="grid md:grid-cols-2 gap-6">
        {useCases.map((useCase, index) => (
          <div
            key={index}
            className="p-8 rounded-2xl bg-gradient-to-br from-white/5 to-white/10 border border-white/10 hover:border-purple-500/50 transition group"
          >
            <div className="flex items-start space-x-4">
              <div className="p-2 bg-purple-500/20 rounded-lg">
                <Check className="w-6 h-6 text-purple-400" />
              </div>
              <div>
                <h3 className="text-xl font-semibold mb-2 group-hover:text-purple-400 transition">
                  {useCase}
                </h3>
                <p className="text-gray-400 text-sm">
                  Accelerate discovery with precise quantum mechanical calculations
                </p>
              </div>
            </div>
          </div>
        ))}
      </div>
    </section>
  );
}