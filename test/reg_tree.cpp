#include <string>
#include <iostream>
#include <iterator>

#include <geomc/Tree.h>
#include <geomc/random/RandomTools.h>

using namespace geom;
using namespace std;


// typedef Tree<const char*, index_t> tree;
// typedef Subtree<const char*, index_t> subt;
// typedef ConstSubtree<const char*, index_t> csubt;


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


void print_subtree(Subtree<const char*, index_t> p, index_t depth) {
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


int main(int argc, char** argv) {
    Tree<const char*, index_t> t;
    
    std::cout << "===== static small tree test =====" << std::endl;
    fill_tree(t.root());
    print_subtree(t.root(), 0);
    std::cout << std::endl;
    
    std::cout << "===== flat tree test =====" << std::endl;
    Tree<const char*, index_t> t2;
    auto r = t2.root();
    *r = "blerg";
    fill_leaf_with_items(r, 30);
    print_subtree(r, 0);
    std::cout << std::endl;
    
    std::cout << "===== simple split() test =====" << std::endl;
    auto c = r.insert_child_node("dingus");
    c.split(compare, 0, "doot");
    print_subtree(r, 0);
    std::cout << std::endl;
    
    std::cout << "===== randomly filled tree =====" << std::endl;
    Tree<const char*, index_t> t3;
    *(t3.root()) = "zoinks";
    populate_tree(t3.root(), 30);
    print_subtree(t3.root(), 0);
    std::cout << std::endl;
    
    std::cout << "===== make a copy of a tree =====" << std::endl;
    Tree<const char*, index_t> t4(t3);
    print_subtree(t4.root(), 0);
    assert(t4 == t3);
    std::cout << std::endl;
    
    return 0;
}

