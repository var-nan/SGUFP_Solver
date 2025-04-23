//
// Created by nandgate on 4/22/2025.
//

#include "../DD.h"
#include <vector>
#include <benchmark/benchmark.h>

using namespace std;

static vector<Inavap::Node> randomNodes(size_t n) {
    vector<Inavap::Node> nodes;
    nodes.reserve(n);
    for (int i = 0; i < n; i++) {
        Inavap::Node node;
        node.globalLayer = n-i;
        nodes.push_back(node);
    }
    return nodes;
}

static void BMVecObject(benchmark::State& state) {
    // create vector of objects.
    for (auto _ :state) {
        state.PauseTiming();
        vector<Inavap::Node> nodes = randomNodes(state.range(0));
        state.ResumeTiming();
        std::sort(nodes.begin(), nodes.end(), [](const Inavap::Node& a, const Inavap::Node& b) {
            return a.globalLayer < b.globalLayer;
        });
    }
}

static void BMVecObjectInsert(benchmark::State& state) {
    for (auto _ :state) {
        state.PauseTiming();
        vector<Inavap::Node> nodes = randomNodes(state.range(0));
        vector<Inavap::Node> actual;
        actual.reserve(nodes.size());
        state.ResumeTiming();
        for (int i = 0; i < nodes.size(); i++) {
            // Inavap::Node node = nodes[i];
            actual.push_back(nodes[i]);
        }
    }
}
static void BMVecObjectPop(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        vector<Inavap::Node> nodes = randomNodes(state.range(0));
        state.ResumeTiming();
        while (!nodes.empty()) {
            Inavap::Node node = nodes.back();
            benchmark::DoNotOptimize(node.globalLayer);
            nodes.pop_back();
        }
    }
}

static vector<Inavap::Node *> randomPointers(size_t n) {
    vector<Inavap::Node *> pointers;
    pointers.reserve(n);
    for (int i = 0; i < n; i++) {
        Inavap::Node *node = new Inavap::Node();
        node->globalLayer = n-i;
        pointers.push_back(node);
    }
    return pointers;
}

static void BMVecPointer(benchmark::State& state) {
    for (auto _ :state) {
        state.PauseTiming();
        vector<Inavap::Node *> pointers = randomPointers(state.range(0));
        state.ResumeTiming();
        std::sort(pointers.begin(), pointers.end(), [](const Inavap::Node *a, const Inavap::Node *b) {
            return a->globalLayer < b->globalLayer;
        });
        state.PauseTiming();
        for (int i = 0; i < pointers.size(); i++) {
            // delete pointers[i];
        }
        state.ResumeTiming();
    }
}

static void BMVecPointerInsert(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        vector<Inavap::Node *> pointers = randomPointers(state.range(0));
        vector<Inavap::Node *> actual;
        actual.reserve(pointers.size());
        state.ResumeTiming();
        for (int i = 0; i < pointers.size(); i++) {
            actual.push_back(pointers[i]);
        }

        state.PauseTiming();
        for (int i = 0; i < pointers.size(); i++) {
            delete pointers[i];
            delete actual[i];
        }
        state.ResumeTiming();
    }
}

static void BMVecPointerPop(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        vector<Inavap::Node *> pointers = randomPointers(state.range(0));
        state.ResumeTiming();

        while (!pointers.empty()) {
            Inavap::Node *node = pointers.back();
            benchmark::DoNotOptimize(node->globalLayer);
            pointers.pop_back();
        }

        state.PauseTiming();
        for (int i = 0; i < pointers.size(); i++) delete pointers[i];
        state.ResumeTiming();
    }
}
BENCHMARK(BMVecObject)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);
BENCHMARK(BMVecPointer)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);

BENCHMARK(BMVecObjectInsert)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);
BENCHMARK(BMVecPointerInsert)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);

BENCHMARK(BMVecObjectPop)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);
BENCHMARK(BMVecPointerPop)->RangeMultiplier(2)->Range(1<<6, 1<<10)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
