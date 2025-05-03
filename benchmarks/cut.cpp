//
// Created by nandgate on 3/25/2025.
//

#include "cut.h"


static void BMCutLinear(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        size_t n = state.range(0);
        double percentage = state.range(1)/10.0;
        auto coeff = generateCut(n);
        auto queries = generateQueries(coeff, percentage);
        Inavap::Cut maincut(200.0, coeff);
        vector<double> answers;
        answers.reserve(queries.size());
        state.ResumeTiming();

        for (auto query : queries) {
            answers.push_back(maincut.get(query));
        }
    }
}

static void BMCutOptimized(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        size_t n = state.range(0);
        double percentage = state.range(1)/10.0;
        auto coeff = generateCut(n);
        auto queries = generateQueries(coeff, percentage);
        NewCut cut (coeff, 200);
        vector<double> answers;
        answers.reserve(queries.size());
        state.ResumeTiming();

        for (auto query : queries) {
            answers.push_back(cut.get(query));
        }

        state.PauseTiming();
        // validate?
        vector<double> actualAnswers;
        Inavap::Cut oldCut(200.0, coeff);
        actualAnswers.reserve(answers.size());
        for (int i = 0; i < queries.size(); i++) {
            double actual = oldCut.get(queries[i]);
            if (actual != answers[i]) {cout << "Wrong answer" << endl; exit(1);}
        }
        state.ResumeTiming();
    }
}

BENCHMARK(BMCutLinear)->Args({1<<6, 5})
        ->Args({1<<7, 5})->Args({1<<8, 5})
        ->Args({1<<9, 5})->Args({1<<10, 5})
        ->Args({1<<6, 9})
        ->Args({1<<7, 9})->Args({1<<8, 9})
        ->Args({1<<9, 9})->Args({1<<10, 9})->Unit(benchmark::kMillisecond);

BENCHMARK(BMCutOptimized)->Args({1<<6, 5})
        ->Args({1<<7, 5})->Args({1<<8, 5})
        ->Args({1<<9, 5})->Args({1<<10, 5})
        ->Args({1<<6, 9})
        ->Args({1<<7, 9})->Args({1<<8, 9})
        ->Args({1<<9, 9})->Args({1<<10, 9})->Unit(benchmark::kMillisecond);

/* based on the benchmark results, BMOptimized cut is significantly faster (atleast 6x in bigger cuts) than
 *  the previous version of the cut */
// BENCHMARK(BMCutLinear)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);
// BENCHMARK(BMCutOptimized)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();