import React, { useState, useCallback } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Card, CardHeader, CardTitle, CardContent } from '@/components/ui/card';

const PoissonCalculator = () => {
  const [params, setParams] = useState({
    copyNumber: 30,
    sampleSize: 8000,
    totalSpecies: 350000,
    meanCopies: 30
  });

  const calculatePoissonSamples = useCallback((lambda, sampleSize, totalSpecies = 350000, meanCopies = 30) => {
    const factorial = (n) => {
      if (n === 0) return 1;
      let result = 1;
      for (let i = 1; i <= n; i++) result *= i;
      return result;
    };

    let sumExp = -Infinity;
    const maxK = Math.min(1000, lambda * 5);
    
    for (let k = 0; k <= maxK; k++) {
      const logPoisson = k * Math.log(lambda) - lambda - Math.log(factorial(k));
      const probMiss = Math.pow(1 - k/(totalSpecies * meanCopies), sampleSize);
      const logProbMiss = Math.log(Math.max(probMiss, 1e-308));
      const logTerm = logPoisson + logProbMiss;
      sumExp = Math.log(Math.exp(sumExp) + Math.exp(logTerm));
    }
    
    return Math.ceil(Math.log(0.05) / sumExp);
  }, []);

  const calculateCaptureProb = useCallback((copyNum, sampleSize, numSamples, totalSpecies, meanCopies) => {
    // For each sample, probability of missing is (1 - k/(N*M))^S
    const probMissSingle = Math.pow(1 - copyNum/(totalSpecies * meanCopies), sampleSize);
    // Probability of missing in all samples
    const probMissAll = Math.pow(probMissSingle, numSamples);
    // Return probability of capturing in at least one sample
    return (1 - probMissAll) * 100;
  }, []);

  const handleCalculate = useCallback((e) => {
    e.preventDefault();
    
    // Calculate required samples
    const samplesNeeded = calculatePoissonSamples(
      params.copyNumber,
      params.sampleSize,
      params.totalSpecies,
      params.meanCopies
    );

    // Generate probability curve data
    const probData = [];
    const maxCopyNum = params.copyNumber * 2;  // Show up to 2x the target copy number
    
    for (let copyNum = 1; copyNum <= maxCopyNum; copyNum++) {
      const prob = calculateCaptureProb(
        copyNum,
        params.sampleSize,
        samplesNeeded,
        params.totalSpecies,
        params.meanCopies
      );
      
      probData.push({
        copyNumber: copyNum,
        probability: prob
      });
    }

    setResult({
      samplesNeeded,
      probData
    });
  }, [params, calculatePoissonSamples, calculateCaptureProb]);

  const [result, setResult] = useState(null);

  const handleInputChange = useCallback((e) => {
    const { name, value } = e.target;
    setParams(prev => ({
      ...prev,
      [name]: Number(value)
    }));
  }, []);

  return (
    <div className="w-full max-w-4xl mx-auto p-4">
      <Card>
        <CardHeader>
          <CardTitle>Poisson Sample Size Calculator</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            <form onSubmit={handleCalculate} className="space-y-6">
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium mb-1">
                    Minimum Copy Number per Species (λ)
                  </label>
                  <input
                    type="number"
                    name="copyNumber"
                    value={params.copyNumber}
                    onChange={handleInputChange}
                    min="1"
                    className="w-full p-2 border rounded"
                    required
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-1">
                    Phage per Sample
                  </label>
                  <input
                    type="number"
                    name="sampleSize"
                    value={params.sampleSize}
                    onChange={handleInputChange}
                    min="1"
                    className="w-full p-2 border rounded"
                    required
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-1">
                    Total Number of Species
                  </label>
                  <input
                    type="number"
                    name="totalSpecies"
                    value={params.totalSpecies}
                    onChange={handleInputChange}
                    min="1"
                    className="w-full p-2 border rounded"
                    required
                  />
                </div>

                <div>
                  <label className="block text-sm font-medium mb-1">
                    Mean Copies per Species
                  </label>
                  <input
                    type="number"
                    name="meanCopies"
                    value={params.meanCopies}
                    onChange={handleInputChange}
                    min="1"
                    className="w-full p-2 border rounded"
                    required
                  />
                </div>
              </div>

              <button
                type="submit"
                className="w-full bg-blue-500 text-white py-2 px-4 rounded hover:bg-blue-600 transition-colors"
              >
                Calculate
              </button>

              {result && (
                <div className="mt-6 p-4 bg-gray-50 rounded-lg">
                  <h3 className="text-lg font-medium mb-2">Results</h3>
                  <p className="text-2xl font-bold text-blue-600">
                    {result.samplesNeeded.toLocaleString()} samples needed
                  </p>
                  <p className="mt-2 text-sm text-gray-600">
                    This will give you a 95% probability of capturing species with {params.copyNumber} or more copies.
                  </p>
                </div>
              )}
            </form>

            {result && (
              <div className="h-96">
                <h3 className="text-lg font-medium mb-4">Detection Probability by Copy Number</h3>
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart
                    data={result.probData}
                    margin={{ top: 20, right: 30, left: 20, bottom: 40 }}
                  >
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis
                      dataKey="copyNumber"
                      label={{
                        value: 'Copy Number',
                        position: 'bottom',
                        offset: 20
                      }}
                    />
                    <YAxis
                      label={{
                        value: 'Detection Probability (%)',
                        angle: -90,
                        position: 'insideLeft',
                        offset: 0
                      }}
                      domain={[0, 100]}
                    />
                    <Tooltip
                      formatter={(value) => [`${value.toFixed(1)}%`, 'Probability']}
                      labelFormatter={(value) => `${value} copies`}
                    />
                    <Line
                      type="monotone"
                      dataKey="probability"
                      stroke="#2563eb"
                      strokeWidth={2}
                      dot={false}
                    />
                    {/* Add reference line at 95% */}
                    <Line
                      type="monotone"
                      dataKey={() => 95}
                      stroke="#dc2626"
                      strokeDasharray="3 3"
                      strokeWidth={1}
                      dot={false}
                    />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            )}
          </div>

          <div className="mt-8 space-y-4 text-sm text-gray-600">
            <div className="bg-gray-50 p-4 rounded-lg">
              <div className="mb-2 font-medium">Formula Used:</div>
              <div className="font-mono text-xs">
                n = ⌈ln(0.05) / ln(Σ[(λᵏe⁻ᵏ/k!) × (1 - k/(N×M))ˢ])⌉
              </div>
              <div className="mt-2 text-xs">
                Probability of detection in n samples = 1 - (1 - k/(N×M))^(S×n)
              </div>
            </div>
            
            <div>
              <p className="font-medium mb-1">Where:</p>
              <ul className="list-disc pl-6 space-y-1">
                <li>λ (lambda): Minimum copy number per species</li>
                <li>N: Total number of species</li>
                <li>M: Mean copies per species</li>
                <li>S: Sample size (phage per sample)</li>
                <li>n: Number of samples needed</li>
                <li>k: Copy number (for probability curve)</li>
              </ul>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
};

export default PoissonCalculator;