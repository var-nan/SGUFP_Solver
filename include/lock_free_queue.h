//
// Created by nandgate on 5/29/2025.
//

#ifndef LOCK_FREE_QUEUE_H
#define LOCK_FREE_QUEUE_H

#include <iostream>
#include <atomic>
#include "../DD.h"
#include <tuple>
#include <cassert>

#define _queue_limit_ 32


typedef Inavap::Node lf_node;
typedef struct {
    lf_node *start;
    lf_node *end;
    size_t n;
} llist;

class master_queue {
    lf_node* head = nullptr;
    size_t size = 0;
public:
    master_queue() = default;

    void push(const llist& nodes) { // bunch of nodes.
        // push to front
        lf_node *start = nodes.start, *end = nodes.end;
        assert((start) && (end));
        size_t n = nodes.n;
        end->next = head;
        head = start;
        size += n;
    }


    /**
     * Returns n nodes from the master's node queue. The return type is a tuple of three
     * elements aliased as llist representing the start and end pointers to 'n' nodes.
     *
     * @param k : number of nodes to pop from the queue.
     * @return llist - tuple(*start, *end, n)
     */
    llist pop_k(size_t k) { // pop 'k' nodes from the front of queue.

        lf_node *curr = head;
        size_t i = 0;

        for (; i < k && curr->next; curr= curr->next) i++;
        llist  result = {head, curr, i};

        head = curr->next;
        curr->next = nullptr;

        size_t j = 0;
        for (const lf_node *start = result.start; start ; start = start->next) j++;
        result.n = j;
        size -= j;
        return result;
    }

    [[nodiscard]] size_t get_size() const noexcept { return size; }

    [[nodiscard]] bool empty() const noexcept { return !size; }
};

class lf_queue {
    std::atomic<size_t> size = 0;
    std::atomic<lf_node *> head = nullptr;

public:
    lf_queue() = default;

    [[nodiscard]] size_t getSize() const noexcept {return size.load(std::memory_order_acquire); }
    [[nodiscard]] bool empty() const noexcept {return size.load(std::memory_order::acquire) == 0; }

    /**
     * Inserts the given nodes to the front of the queue.
     * The function can be called by both the master and the worker.
     */
    void push(const llist &nodes) {

        /* At any point of time, only one thread (either master or worker) will call this function.
         * When the master is calling this function, the worker must be in waiting state */
        lf_node *start = nodes.start, *end = nodes.end;
        size_t n = nodes.n;
        end->next = head.load(std::memory_order::relaxed);
        head = start;
        head.store(start, std::memory_order::release); // could be relaxed.
        size.fetch_add(n, std::memory_order::release); // cannot convince this memory order.
    }

    /**
     * Pops the top most node from the queue and returns its pointer. Effectively decrements
     * the queue's size by 1.
     */
    lf_node *pop() {

        if (!head) return nullptr;
        lf_node *rv = head.load(std::memory_order::acquire);// can be relaxed.
        head.store(rv->next); // change to release or relaxed later.
        size.fetch_sub(1, std::memory_order::acq_rel);
        rv->next = nullptr;
        return rv;
    }

    /**
     * The function returns nodes from the worker queue.
     * Should be called by the master.
     * @param proportion - proportion of nodes to pop from worker queue.
     * @return llist {start,end, n}
     */
    [[gnu::optimize("O0")]] llist m_pop(double proportion) {

        proportion = 1 - proportion; // % of nodes left in queue after successful completion of this operation.

        size_t sz = size.load(); // memory_order::acquire
        if (sz < _queue_limit_) return {nullptr, nullptr, 0};

        size_t n_skip = static_cast<size_t>(static_cast<double>(sz) * proportion);// number of nodes to skip from head.
        size_t rem = sz - n_skip; // number of nodes popped from the queue from the end.
        size_t k = n_skip;

        lf_node *start = head; // need SC here.

        while (--n_skip && start->next) start = start->next;

        size_t ssz = size.load(); // memory order: acquire
        if (ssz <= (sz - (k>>1))) {
            return {nullptr, nullptr, 0}; // nodes are popped very quick, retry.
        }

        lf_node *begin = start->next;
        start->next = nullptr;
        size.fetch_sub(rem);

        // reach last of the list.
        lf_node *end = begin;
        while (end->next) end = end->next;

        return {begin, end,rem};
    }
};

#endif //LOCK_FREE_QUEUE_H
