#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Tree

#include <boost/test/unit_test.hpp>
#include <geomc/Tree.h>
#include <geomc/random/MTRand.h>
#include <geomc/CircularBuffer.h>
#include <geomc/shape/Rect.h>

// todo:
// check different sizes of tree
// multiple random queries
// queries that are likely to succeed
// queries that are likely to fail
// check that item count actually matches item count
//   same with node count
//   (wrap that into an integrity_check())
//   rm tree-printing
// try removing random nodes
// try emptying a tree by removing random nodes one at a time.
// try with heap memory objects (e.g. behind a shared_ptr)


using namespace geom;
using namespace std;

#define RANDOM_SEED 493866457312

MTRand rng = MTRand(RANDOM_SEED);


inline index_t compare(index_t a, index_t b) {
    return a - b;
}


void fill_tree(const Subtree<const char*, index_t>& t) {
    *t = "i am root";
    // three kids for root
    auto c1  = t.insert_child_node("c1");
    auto c2  = t.insert_child_node("c2");
    auto c3  = t.insert_child_node("c3");
    // one grandchild
    auto gc1 = c2.insert_child_node("gc1");
    // items. items for you.
    c1.insert_item(1);
    c1.insert_item(2);
    c1.insert_item(3);
    gc1.insert_item(4);
    gc1.insert_item(5);
    gc1.insert_item(6);
}


template <typename T>
void fill_leaf_with_items(Subtree<T, index_t> subtree, index_t n) {
    for (index_t i = 0; i < n; ++i) {
        auto b = subtree.items_begin();
        if (subtree.item_count() > 0) {
            index_t j = rng.rand(subtree.item_count()); // xxx need to be able to insert at end
            std::advance(b, j);
        }
        subtree.insert_item(b, rng.rand<index_t>());
    }
}


char* random_name() {
    unsigned char i = rng.rand<int>() & 0xff;
    std::string name = std::to_string(i);
    size_t l = name.length();
    char* s = new char[l + 1];
    std::copy(name.begin(), name.end(), s);
    s[l] = '\0';
    return s;
}


void randomly_place_item(const Subtree<const char*, index_t>& tree) {
    if (tree.node_count() > 0) {
        // descend into a random child
        auto child = tree.begin();
        std::advance(child, rng.rand(tree.node_count()));
        randomly_place_item(child);
    } else if (tree.item_count() > 8) {
        index_t piv  = (*tree.items_end() + *tree.items_begin()) / 2;
        auto left_c  = tree.insert_child_node(random_name());
        auto right_c = left_c.split(compare, piv, random_name());
        // try again
        randomly_place_item(tree);
    } else {
        // place it here.
        tree.insert_item(rng.rand<index_t>());
    }
}


void populate_tree(const Subtree<const char*, index_t>& tree, index_t n) {
    for (index_t i = 0; i < n; ++i) {
        randomly_place_item(tree);
    }
}


template <typename T, index_t N>
std::ostream& operator<<(std::ostream& s, const Rect<T,N>& r) {
    s << "[" << r.min() << " : " << r.max() << "]";
    return s;
} 


template <typename T>
void print_subtree(Subtree<T, index_t> p, index_t depth=0) {
    cout << std::string(depth * 2, ' ') << "(+) " << *p << " ";
    if (p.node_count() == 0) {
        cout << p.item_count() << endl;
        auto indent = std::string(2 * (depth + 1), ' ');
        for (auto i = p.items_begin(); i != p.items_end(); ++i) {
            cout << indent << *i << endl;
        }
    } else {
        cout << p.node_count() << endl;
        for (auto c = p.begin(); c != p.end(); ++c) {
            print_subtree(c, depth + 1);
        }
    }
}

void split_subtree(Subtree<Rect<index_t, 1>, index_t> n) {
    const index_t max_arity = 4;
    if (n.item_count() > max_arity) {
        // make a new child identical to ourselves; put it on the stack
        CircularBuffer<Subtree<Rect<index_t, 1>, index_t>, 8> buf;
        buf.push_back(n.insert_child_node(*n));
        
        // split into siblings until we are full, or all the siblings are small enough.
        while (buf.size() > 0 and n.node_count() < max_arity) {
            auto child = buf.pop_front();
            if (child.item_count() > max_arity) {
                double avg = 0;
                index_t ct = 0;
                for (auto k = child.items_begin(); k != child.items_end(); ++k) {
                    avg += *k;
                    ++ct;
                }
                avg /= ct;
                auto new_node = child.split(compare, avg);
                if (new_node == child.end()) {
                    continue;
                }
                // if the split produced an empty node,
                // don't try to split again to avoid an infinite loop.
                if (new_node.item_count() == 0) {
                    n.erase(new_node);
                } else if (child.item_count() == 0) {
                    n.erase(child);
                } else {
                    // both nodes are potential candidates for further splitting.
                    buf.push_back(new_node);
                    buf.push_back(child);
                }
            }
        }
        
        if (n.node_count() > 1) {
            // subdivide / recurse on new children
            for (auto child = n.begin(); child != n.end(); ++child) {
                split_subtree(child);
            }
        } else {
            // split was unsuccessful.
            n.erase(n.begin());
        }
    }
    
    // rebound this node
    // not optimal but meh
    Rect<index_t, 1> bnd;
    for (auto i = n.items_begin(); i != n.items_end(); ++i) {
        bnd |= *i;
    }
    *n = bnd;
}


bool bound(Subtree<Rect<index_t, 1>, index_t>& t, index_t s) {
    return t->contains(s);
}


void populate_search_tree(Tree<Rect<index_t, 1>, index_t>& tree, index_t n) {
    auto r = tree.root();
    *r = Rect<index_t, 1>();
    fill_leaf_with_items(r, n);
    // split ourselves
    split_subtree(r);
    // test a buncha queries
    for (index_t j = 0; j < 200; ++j) {
        const index_t query_pt = rng.rand<index_t>();
        index_t ct = 0;
        for (auto i = r.query(query_pt, bound, bound); i != r.end() and ct < n; ++i, ++ct) {
            auto bnd = **i;
            BOOST_CHECK(bnd.contains(query_pt));
        }
    }
}


///////////// test suites /////////////


BOOST_AUTO_TEST_SUITE(tree_construction)


BOOST_AUTO_TEST_CASE(construct_small_tree) {
    Tree<const char*, index_t> t;
    fill_tree(t.root());
    BOOST_CHECK(t.size() == 5);
    BOOST_CHECK(t.item_count() == 6);
}


BOOST_AUTO_TEST_CASE(split_flat_tree) {
    const index_t count = 30;
    Tree<const char*, index_t> t;
    auto r = t.root();
    *r = "blerg";
    fill_leaf_with_items(r, count);
    
    // verify all items are accessible.
    auto i = r.items_begin();
    for (index_t j = 0; j < count; ++j, ++i) {
        index_t b = *i;
        BOOST_REQUIRE(b == b);
    }
    
    // add a single child:
    //   r
    //   |
    //   c
    
    auto c = r.insert_child_node("dingus");
    // does the new child inherit the parent's items?
    BOOST_CHECK(c.item_count() == count);
    // is the parent of the new child the root?
    BOOST_REQUIRE(c.parent() == r);
    // is the parent's first child the new node?
    BOOST_REQUIRE(r.begin() == c);
    
    // now test splitting:
    //    r
    //   / \
    //  s   c
    auto s = c.split(compare, 0, "doot");
    // parent has one more child?
    BOOST_REQUIRE(r.node_count() == 2);
    // is the new node's value as specified?
    BOOST_CHECK(strcmp(*s, "doot") == 0);
    // is the new node still a child of its parent?
    BOOST_REQUIRE(s.parent() == r);
    // the new node should have been inserted before `c`, the first child.
    BOOST_REQUIRE(r.begin() == s);
    // split nodes are siblings:
    BOOST_REQUIRE(std::next(s) == c);
    // item count makes sense:
    BOOST_CHECK(r.item_count() == s.item_count() + c.item_count());
    // items are contiguous:
    BOOST_CHECK(s.items_end() == c.items_begin());
    // parent's edge items are correct:
    BOOST_REQUIRE(r.items_begin() == s.items_begin());
    BOOST_REQUIRE(r.items_end() == c.items_end());
    
    print_subtree(r);
}


BOOST_AUTO_TEST_CASE(copy_deep_tree) {
    Tree<const char*, index_t> t;
    *(t.root()) = "zoinks";
    populate_tree(t.root(), 30);
    print_subtree(t.root());
    
    Tree<const char*, index_t> t2(t);
    BOOST_CHECK(t2 == t);
}


BOOST_AUTO_TEST_CASE(search_tree) {
    Tree<Rect<index_t, 1>, index_t> t;
    // fuzz this baby
    for (index_t i = 0; i < 1024; ++i) {
        populate_search_tree(t, rng.rand<index_t>(200));
        t.root().clear();
        BOOST_CHECK_EQUAL(t.size(), 1);
        BOOST_CHECK_EQUAL(t.item_count(), 0);
    }
}


BOOST_AUTO_TEST_SUITE_END()
