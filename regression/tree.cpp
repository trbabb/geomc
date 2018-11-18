#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Hello

#include <boost/test/unit_test.hpp>
#include <geomc/Tree.h>
#include <geomc/random/RandomTools.h>

using namespace geom;
using namespace std;


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


void fill_leaf_with_items(Subtree<const char*, index_t> tree, index_t n) {
    Random* rng = geom::getRandom();
    for (index_t i = 0; i < n; ++i) {
        auto b = tree.items_begin();
        if (tree.item_count() > 0) {
            std::advance(b, rng->rand(tree.item_count()));
        }
        tree.insert_item(b, rng->rand<index_t>());
    }
}


char* random_name(Random* rng) {
    unsigned char i = rng->rand<int>() & 0xff;
    std::string name = std::to_string(i);
    size_t l = name.length();
    char* s = new char[l + 1];
    std::copy(name.begin(), name.end(), s);
    s[l] = '\0';
    return s;
}


void randomly_place_item(const Subtree<const char*, index_t>& tree) {
    Random* rng = geom::getRandom();
    if (tree.node_count() > 0) {
        // descend into a random child
        auto child = tree.begin();
        std::advance(child, rng->rand(tree.node_count()));
        randomly_place_item(child);
    } else if (tree.item_count() > 8) {
        index_t piv  = (*tree.items_end() + *tree.items_begin()) / 2;
        auto left_c  = tree.insert_child_node(random_name(rng));
        auto right_c = left_c.split(compare, piv, random_name(rng));
        // try again
        randomly_place_item(tree);
    } else {
        // place it here.
        tree.insert_item(rng->rand<index_t>());
    }
}


void populate_tree(const Subtree<const char*, index_t>& tree, index_t n) {
    for (index_t i = 0; i < n; ++i) {
        randomly_place_item(tree);
    }
}


void print_subtree(Subtree<const char*, index_t> p, index_t depth=0) {
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

BOOST_AUTO_TEST_SUITE(tree_construction)

BOOST_AUTO_TEST_CASE(construct_small_tree)
{
    Tree<const char*, index_t> t;
    fill_tree(t.root());
    BOOST_CHECK(t.size() == 5);
    BOOST_CHECK(t.item_count() == 6);
}

BOOST_AUTO_TEST_CASE(split_flat_tree)
{
    const index_t count = 30;
    Tree<const char*, index_t> t;
    auto r = t.root();
    *r = "blerg";
    fill_leaf_with_items(r, count);
    
    // verify all items are accessible.
    auto i = r.items_begin();
    for (index_t j = 0; j < count; ++j, ++i) {
        index_t b = *i;
        BOOST_REQUIRE(i == i);
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

BOOST_AUTO_TEST_CASE(copy_deep_tree)
{
    Tree<const char*, index_t> t;
    *(t.root()) = "zoinks";
    populate_tree(t.root(), 30);
    print_subtree(t.root());
    
    Tree<const char*, index_t> t2(t);
    BOOST_CHECK(t2 == t);
}

BOOST_AUTO_TEST_SUITE_END()
