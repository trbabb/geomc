#include <boost/iterator/iterator_facade.hpp>
#include <geomc/geomc_defs.h>
#include <geomc/Templates.h>
#include <list>
#include <type_traits>

// todo: t.query<Key, BoundingFn(k,node), MatchingFn(k,item)>(Key k) -> item_iterator (derived)
//       - or t.visit(visitor)
//         - visitor has visit(node)
//         - visitor has should_traverse(node)
//         > the former (iterator) is more flexible with less infrastructure.
//           the traversal is interchangeable and generic. with a visitor you
//           have to re-implement the traversal, and also bundle your logic 
//           into a class or a functor, instead of doing it in-place. too fussy.
//       - something similar for trees with no items.
//         - t.query<Key, BoundingFn(k,node), MatchingFn(k,node)>(Key k) -> node_iterator (derived)
//       *** consider what happens if the tree is mutated underneath you
// todo: permit key/value structure. use BoundingFn for search.
// todo: permit a distance metric; allow KNN.
//       >> implement as a free function.
// todo: consider factoring the above simply as free functions which implement algorithms.
// todo: specialize for when leaf items or node items are void
// todo: smarter/more optimized underlying data structure
// todo: "flat array" variant tree. abstract NodeRef and ItemRef to be ptrs.
// todo: might be good to implement take(subtree) which removes 
//       the subtree from its source and adopts it. if it belongs to the
//       same tree, a copy could be avoided. be sure to permit the case
//       of re-rooting (i.e. root().take(some_descendent_of_self)).
//       generally only cheap if `list` permits slicing and splicing.
// todo: it is not possible to have an empty tree. (can have zero items, but must have one node).

// todo: get the docs correct about what end() and begin() iterators are invalidated.
//       -> begin() on this node if we insert before begin
//       -> items_begin() on this node and any parent if we insert before begin
//       -> {items,}end() on predecessors / their parents if we insert before begin
//       -> {items,}begin() and {items,}end() on this node if it's the first insertion
//       -> more for deletion, splitting, adopting, (taking)
//       all this is complicated, so editing while using a range for() is
//       just all around a bad idea.

// xxx: how to subclass?
//   > subclass Tree and call it, e.g. AssociativeTree and/or RTree
//   > subclass the two iterators
//     > may have to templatize SubtreeBase over Tree type so that member and
//       typedefs are in agreement.
//     > SubtreeBase also directly refers to Subtree and ConstSubtree, which
//       will have to be corrected somehow.
//   > have the subclassed Tree return the new iterators from root() and begin()/end()
//   - might need a way to chain delegation of self_t up the iterator inheritance hierarchy.
//     - tree defines what its iterators are, iterators are passed tree_t, iterators
//       ask tree_t what they are?
//   > or! the overridden iterator just subclasses from SubtreeBase<...>
//     since the mutators are going to be different.
//     - and to mutate, we internally construct a Subtree<> from _root,
//       but never become one or hand it out?
//   - note with this scheme algorithms have to pick the subclass of tree
//     they work on :|
//     > ...which is probably OK.

// ...OR:
//   - an R-tree "has a" Tree, not "is a" tree.
//   - R-tree has insert(item), erase({item,node}_iterator), find(item), find_parent(item)
//     which return iterators, and call the tree iterators.
//     - uses Bound(item) -> Boundary
//     - uses Bound(Boundary, item) -> bool
//     - uses Union(Boundary, Boundary) -> Boundary
//   - possibly only let the user have a const_subtree iterator into the master tree.
//   - consider how to make ListTree and ArrayTree interchangeable within

// todo: formally define the sorting order


namespace geom {

    
/** @addtogroup storage
 *  @{
 */


// forward decls
template <typename NodeItem, typename LeafItem, bool Const> class SubtreeBase;
template <typename NodeItem, typename LeafItem> class Subtree;
template <typename NodeItem, typename LeafItem> class ConstSubtree;


/**
@brief A dynamic tree of arbitrary arity.

Each node in the tree may store an associated `NodeItem` object. 
A node may have any number of child nodes.

Leaf nodes may own a list of `LeafItem` objects. `LeafItem`s are kept
in (contiguous) depth-first order, and internal nodes may therefore provide 
iterators to the first and last `LeafItem`s in their respective subtrees.

All mutation and access to a tree is provided through the Subtree and 
ConstSubtree classes; both thin iterators pointing into the Tree data structure.

Trees are kept in "sibling-contiguous" order, as replicated by this pseudocode:
 
<pre>
traverse(node):
    for each child c of node:
        visit(c)
    for each child c of node:
        traverse(c)
</pre>

which produces nodes in this order:

<pre>
[node] ... [c1 c2 c3 c4 ...]
    [children of c1]
        [descendents of c1...]
    [children of c2]
        [descendents of c2...]
    ...
</pre>

In this manner, a subtree may be distant from its root, but that entire
subtree is contiguous. This also preserves the contiguousness of siblings.

@tparam NodeItem Type of data to be kept with internal tree nodes.
@tparam LeafItem Type of data to be kept by leaf nodes.
*/
template <typename NodeItem, typename LeafItem>
class Tree {
    
    friend class SubtreeBase<NodeItem, LeafItem, true>;
    friend class SubtreeBase<NodeItem, LeafItem, false>;
    friend class ConstSubtree<NodeItem, LeafItem>;
    friend class Subtree<NodeItem, LeafItem>;
    
    struct Node;
    
    typedef typename std::list<Node>::iterator           NodeRef;
    typedef typename std::list<LeafItem>::iterator       ItemRef;
    typedef typename std::list<LeafItem>::const_iterator ConstItemRef;
    
    struct Node {
        
        Node() {}
        
        template <typename ... Args>
        Node(Tree<NodeItem, LeafItem>* storage,
             NodeRef parent,
             Args&& ... args):
                parent(parent),
                child_first (storage->_nodes.end()),
                child_last  (storage->_nodes.end()),
                items_first (storage->_items.end()),
                items_last  (storage->_items.end()),
                n_items     (0),
                n_children  (0),
                data        (args...) {}
        
        NodeRef  parent;
        NodeRef  child_first;
        NodeRef  child_last;
        ItemRef  items_first;
        ItemRef  items_last;
        
        index_t  n_items;    // i.e. all descendents
        index_t  n_children; // i.e. direct children
        
        NodeItem data;
        
    };
    
    std::list<Node>     _nodes;
    std::list<LeafItem> _items;
    
public:
    
    /// Construct an empty Tree.
    Tree() {
        _nodes.push_back(Node(this, _nodes.end()));
    }
    
    
    /// Default move constructor.
    Tree(Tree<NodeItem, LeafItem>&& other) = default;
    
    
    /**
     * @brief Construct a tree having a single root node, and copy the
     * items in `[begin, end)` to it.
     */
    template <typename LeafItemIterator>
    Tree(LeafItemIterator begin, LeafItemIterator end) {
        // load items
        _items.assign(begin, end);
        
        // create root
        _nodes.push_back(Node(this, _nodes.end()));
        NodeRef node = _nodes.begin();
        
        // assign items to root
        node->items_first = _items.begin();
        node->items_last  = _items.end(); --(node->items_last);
        node->n_items     = _nodes.size();
    }
    
    
    /// Make a copy of the given tree.
    Tree(const Tree<NodeItem, LeafItem>& other):
            _nodes(other._nodes),
            _items(other._items) {
        recalculate_tree();
    }
    
    
    /// Construct a new tree from the given subtree.
    Tree(const ConstSubtree<NodeItem, LeafItem>& other):
            _nodes(other.subtree_begin(), other.subtree_end()),
            _items(other.items_begin(), other.items_end()) {
        // make a root (we only added the children)
        _nodes.push_front(*other._root);
        recalculate_tree();
    }
    
    
    // xxx: Todo: make a constructor which populates the NodeItem
    
    
    /// Get a subtree iterator to the root of this tree.
    inline Subtree<NodeItem, LeafItem> root() {
        return Subtree<NodeItem, LeafItem>(_nodes.begin(), this);
    }
    
    
    /// Get a const subtree iterator to the root of this tree.
    inline ConstSubtree<NodeItem, LeafItem> root() const {
        return ConstSubtree<NodeItem, LeafItem>(_nodes.begin(), this);
    }
    
    
    /** 
     * @brief Return an iterator to the first item in the tree (root), 
     * in sibling-contiguous order.
     */
    Subtree<NodeItem, LeafItem> begin() {
        return root();
    }
    
    
    /**
     * @brief Return an iterator to the last (off-end) item in the tree, 
     * in sibling-contiguous order.
     */
    Subtree<NodeItem, LeafItem> end() {
        return Subtree<NodeItem, LeafItem>(_nodes.end(), this);
    }
    
    
    /**
     * @brief Return a const iterator to the first item in the tree (root), 
     * in sibling-contiguous order.
     */
    ConstSubtree<NodeItem, LeafItem> begin() const {
        return root();
    }
    
    
    /**
     * @brief Return a const iterator to the last (off-end) item in the tree, 
     * in sibling-contiguous order.
     */
    ConstSubtree<NodeItem, LeafItem> end() const {
        return ConstSubtree<NodeItem, LeafItem>(_nodes.end(), this);
    }
    
    
    /**
     * @brief Convert a ConstSubtree into a Subtree by removing its constness.
     *
     * Allows code with mutable access to the source tree to use ConstSubtree
     * handles for any mutation with constant overhead.
     *
     * If the supplied ConstSubtree does not belong to this Tree, then `end()` is returned.
     */
    inline Subtree<NodeItem, LeafItem> subtree(ConstSubtree<NodeItem, LeafItem>& i) {
        if (i._storage == this) {
            return Subtree<NodeItem, LeafItem>(i._root, this);
        } else {
            return end();
        }
    }
    
    
    /// Total number of nodes in this Tree.
    inline size_t size() const {
        // may be slow for dumb implementations of std::list :(
        return _nodes.size();
    }
    
    
    /// Total number of items in this Tree.
    inline size_t item_count() const {
        return _nodes.begin()->n_items;
    }
    
    
    /// Tree equality
    bool operator==(const Tree<NodeItem, LeafItem>& other) const {
        if (this == &other) return true;
        if (item_count() != other.item_count()) return false;
        if (_items       != other._items)       return false;
        
        // because the references will point to different memory,
        // we can't compare raw data. we have to examine the structure of the tree.
        auto a_n =       _nodes.begin();
        auto b_n = other._nodes.begin();
        for (; a_n != _nodes.end() and b_n != other._nodes.end(); ++a_n, ++b_n) {
            if (a_n->n_items    != b_n->n_items)    return false;
            if (a_n->n_children != b_n->n_children) return false;
            if (a_n->data       != b_n->data)       return false;
        }
        
        // must have same total number of nodes; not just a matching prefix subtree.
        // we don't check item_count() beforehand because it may not be stored, 
        // and will thus doubly incur the linear cost we have pay above.
        if (not (a_n == _nodes.end() and b_n == other._nodes.end())) return false;
        
        return true;
    }
    
    
    /// Tree inequality
    bool operator!=(const Tree<NodeItem, LeafItem>& other) const {
        return not (*this == other);
    }
    
    
protected:
    
    inline void recalculate_tree() {
        NodeRef n = _nodes.begin();
        n->child_first = std::next(n);
        n->child_last  = _nodes.end(); --(n->child_last);
        n->items_first = _items.begin();
        n->items_last  = _items.end(); --(n->items_last);
        root().recalculate_references();
    }
    
    
    // called when a leaf node adds or removes one or more of its edge items.
    // the edge items for the anscestors potentially need to be updated to
    // include/exclude the affected items.
    void update_boundary_items(
            NodeRef leaf,
            ItemRef new_first,
            ItemRef new_last,
            index_t delta_items) {
        
        NodeRef prev_child = _nodes.end();
        NodeRef n = leaf;
        
        ItemRef old_first = n->items_first;
        ItemRef old_last  = n->items_last;
        
        // update ancestors' item count and boundary items
        while (n != _nodes.end()) {
            n->n_items += delta_items;
            
            if (n->n_items == 0) {
                n->items_first = n->items_last = _items.end();
            } else {
                // If we ascended from an edge node, that node certainly
                // controls (or used to control) our edge item. 
                // But it might have also controlled our edge item if the 
                // edge sibling doesn't have any items. In this case,
                // its previous edge item will also be ours.
                if (n->child_first == prev_child or 
                        n->items_first == old_first) {
                    n->items_first = new_first;
                }
                if (n->child_last == prev_child or
                        n->items_last == old_last) {
                    n->items_last = new_last;
                }
            }
            
            prev_child = n;
            n = n->parent;
        }
    }
    
};



/**
 * @brief Base class for all iterators into Trees.
 * 
 * A Subtree's lifetime is valid no longer than the Tree from which it came.
 *
 * Mutations to the source Tree generally invalidate end() iterators.
 *
 * It is invalid to dereference, mutate, or increment an end() iterator.
 */
template <typename NodeItem, typename LeafItem, bool Const=true>
class SubtreeBase {
    
protected:
    
    typedef Tree<NodeItem, LeafItem>  Storage;
    typedef typename Storage::NodeRef NodeRef;
    typedef typename Storage::ItemRef ItemRef;
    typedef typename Storage::ItemRef ConstItemRef;
    typedef typename Storage::Node    Node;
    typedef typename std::conditional<
            Const,
            const Storage*,
            Storage*>::type           StorageRef;
    
    
    StorageRef _storage;
    NodeRef    _root;
    
    
    // construct an iterator to a particular node
    SubtreeBase(const NodeRef& root, StorageRef storage):
        _storage(storage),
        _root(root) {}
    
public:
    
    /// Self type. A ConstSubtree if this is a const iterator; a Subtree otherwise.
    typedef typename std::conditional<
            Const,
            ConstSubtree<NodeItem, LeafItem>,
            Subtree<NodeItem, LeafItem>
        >::type self_t;
    /// Const iterator to `LeafItem`s.
    typedef ConstItemRef    const_item_iterator;
    /// A (possibly const) iterator to `LeafItem`s.
    typedef typename std::conditional<
            Const,
            ConstItemRef,
            ItemRef
        >::type item_iterator;
    
    /// A (possibly const) iterator over subtrees.
    typedef self_t iterator;
    /// Const iterator over subtrees
    typedef ConstSubtree<NodeItem, LeafItem>  const_iterator;
    /// The tree's `NodeItem` type
    typedef NodeItem value_type;
    /// A (possibly const) reference to the tree's `NodeItem` type
    typedef typename std::conditional<
            Const,
            const NodeItem&,
            NodeItem&
        >::type reference;
    /// A const reference to the tree's `NodeItem` type
    typedef const NodeItem& const_reference;
    /// A (possibly const) pointer to the tree's `NodeItem` type
    typedef typename std::conditional<
            Const,
            const NodeItem*,
            NodeItem*
        >::type pointer;
    /// Type of tree into which this iterator points.
    typedef Tree<NodeItem, LeafItem> tree_t;
    /// Iterator difference type.
    typedef typename std::iterator_traits<NodeRef>::difference_type difference_type;
    /// Iterator category.
    typedef typename std::iterator_traits<NodeRef>::iterator_category iterator_category;
    
    
    /**
     * @brief A forward iterator over this type of Subtree, which visits 
     * all the nodes which pass the supplied `TestFn`.
     *
     * @tparam I The type of Subtree to visit.
     * @tparam BoundingFn A callable object which accepts an `I` and
     * returns `true` if that Subtree may contain items which pass `TestFn`.
     * @tparam TestFn A callable object which accepts an `I` and returns
     * `true` if that Subtree should be visited by this `QueryIterator`.
     */
    template <typename I, typename BoundingFn, typename TestFn>
    class QueryIterator : public boost::iterator_facade<
            QueryIterator<I, BoundingFn, TestFn>, // curiously recurring!
            I,                                    // type referred to
            boost::forward_traversal_tag> {       // forward-only iterator (for now)
        
        friend class boost::iterator_core_access;
        
        I _item;
        I _end;
        BoundingFn _bound;
        TestFn     _test;
        
    public:
        
        /// Construct a new QueryIterator.
        QueryIterator(const I& start, BoundingFn bound, TestFn test):
                _item(start),
                _end(_item.subtree_end()), 
                _bound(bound),
                _test(test) {
            if (not _bound(_item)) {
                _item = _end;
            } else if (not _test(_item)) {
                increment();
            }
        }
        
        /// Return true if this iterator points at the subtree at `other`.
        bool operator==(const I& other) const {
            return other == _item;
        }
    
    private:
        
        void increment() {
            ++_item;
            while (_item != _end and not _test(_item)) {
                if (_bound(_item)) {
                    // descend into tree
                    _item = _item.begin();
                } else {
                    // jump over tree
                    _item = _item.subtree_end();
                }
            }
        }
        
        bool equal(const QueryIterator& other) const {
            return _item == other._item;
        }
        
        I& dereference() const {
            return _item;
        }
        
    };
    
    /// Convenience test function for `query()` which returns `true` on all nodes.
    static bool visit_all_nodes(const self_t& s)  { return true; }
    /// Convenience test function for `query()` which returns `true` on all nodes with zero child nodes.
    static bool visit_all_leaf_nodes(const self_t& s) { return s.node_count() == 0; }
    
    
    /************************************
     * Methods                          *
     ************************************/

public:
    
    
    /// Obtain the tree into which this iterator points.
    inline typename std::conditional<Const, const tree_t&, tree_t&>::type
    tree() const {
        return _storage;
    }
    
    /**
     * @brief Get the parent node.
     *
     * The parent of `tree()->root()` is `tree()->end()`.
     */
    inline self_t parent() const {
        return self_t(_root->parent, _storage);
    }
    
    /// Return whether this node is the root of the entire tree.
    inline bool is_root() const {
        return _root == _storage->_nodes.begin();
    }
    
    /// `+i`: Become first child
    inline self_t& operator+() {
        _root = _root->child_first;
        return *static_cast<self_t*>(this);
    }
    
    /// `-i`: Become parent
    inline self_t& operator-() {
        _root = _root->parent;
        return *static_cast<self_t*>(this);
    }
    
    /// `++i`: Become next sibling
    inline self_t& operator++() {
        _root = ++_root;
        return *static_cast<self_t*>(this);
    }
    
    /// `--i`: Become previous sibling
    inline self_t& operator--() {
        _root = --_root;
        return *static_cast<self_t*>(this);
    }
    
    /// `i++`: Become next sibling
    inline self_t operator++(int) {
        self_t tmp = *this;
        ++(*this);
        return tmp;
    }
    
    /// `i--`: Become previous sibling
    inline self_t operator--(int) {
        self_t tmp = *this;
        --(*this);
        return tmp;
    }
    
    /// Returns `true` iff `other` points to the same node of the same tree.
    inline bool operator==(const self_t& other) const {
        return _root == other._root;
    }
    
    /// Returns `true` iff `other` does not point to the same node of the same tree.
    inline bool operator!=(const self_t& other) const {
        return _root != other._root;
    }
    
    /// `*i`: Get the value of the current node
    inline reference operator*() const {
        return _root->data;
    }
    
    /// `i->...`: Access member of current node
    inline pointer operator->() const {
        return &_root->data;
    }
    
    /// Number of direct child nodes
    inline index_t node_count() const {
        return _root->n_children;
    }
    
    /// Number of leaf items in this subtree
    inline index_t item_count() const {
        return _root->n_items;
    }
    
    /// Get first child
    inline self_t begin() const {
        return self_t(_root->child_first, _storage);
    }
    
    /// Get last (off-end) child
    inline self_t end() const {
        NodeRef tmp = _root->child_last;
        // cannot increment the ::end() iterator
        return self_t((_root->n_children > 0) ? ++tmp : tmp, _storage);
    }
    
    /// Get first object inside this subtree
    inline item_iterator items_begin() const {
        return _root->items_first;
    }
    
    /**
     * @brief Get last (off-end) object in this subtree. 
     * It is invalid to increment or dereference this iterator.
     */
    inline item_iterator items_end() const {
        item_iterator tmp = _root->items_last;
        // cannot increment the ::end() iterator
        return (_root->n_items > 0) ? ++tmp : tmp;
    }
    
    
    /**
     * @brief Return the first subtree in a sequence covering all the
     * nodes in this subtree, beginning with the first child of this node.
     */
    inline self_t subtree_begin() const {
        return _root->child_first;
    }
    
    
    /**
     * @brief Return the last (off-end) subtree in the sequence of all subtrees
     * below this node.
     */
    inline self_t subtree_end() const {
        return this->node_successor(_root);
    }
    
    
    /**
     * @brief Return an iterator to the first item in the entire tree
     * (the global root).
     *
     * Alias for `this->tree().begin()`.
     */
    inline self_t global_begin() const {
        return self_t(_storage->_nodes.begin(), _storage);
    }
    
    /**
     * @brief Return an iterator to the last item in the entire tree
     * (the global off-end node).
     *
     * Alias for `this->tree().end()`.
     */
    inline self_t global_end() const {
        return self_t(_storage->_nodes.end(), _storage);
    }
    
    
    /**
     * @brief Exhaustively search for the direct parent node of item `i` in 
     * this subtree.
     *
     * If `i` is not in the subtree under this node, then return `end()`.
     */
    self_t find_parent(const const_item_iterator& i) const {
        NodeRef subtree_last = node_successor(_root);
        for (NodeRef n = _root; n != subtree_last; ++n) {
            if (n->n_children == 0 and n->n_items > 0) {
                ItemRef end_item = n->items_last; ++end_item;
                for (ItemRef j = n->items_first; j != end_item; ++j) {
                    if (j == i) return self_t(n, _storage);
                }
            }
        }
        return self_t(_storage->_nodes.end(), _storage);
    }
    
    
    /**
     * @brief Search for the direct parent node of item `i` in 
     * this subtree, using `BoundingFn` to exclude subtrees which 
     * cannot contain `i`.
     *
     * If `i` is not in the subtree under this node, then return `end()`.
     *
     * Because `BoundingFn` is templated, it may be inlined for performance.
     *
     * @param i The item to find the parent of.
     * @tparam BoundingFn A function which returns `true` if node `n` might possibly 
     * contain item `i`.
     */
    template <bool BoundingFn(const NodeItem&, const LeafItem&)>
    self_t find_parent(const const_item_iterator& i) const {
        NodeRef cur_root     = _root;
        NodeRef subtree_last = node_successor(_root);
        while (cur_root != subtree_last) {
            if (BoundingFn(cur_root->data, *i)) {
                // cur_root may contain `i`
                if (cur_root->n_children > 0) {
                    // recurse
                    cur_root = cur_root->child_first;
                } else {
                    // look for item
                    ItemRef end_item = cur_root->items_last; ++end_item;
                    for (ItemRef j = cur_root->items_first; j != end_item; ++j) {
                        // found it
                        if (j == i) return self_t(cur_root, _storage);
                    }
                    // our current node is a leaf, so there is no subtree to jump over.
                    ++cur_root;
                }
            } else {
                // check next sibling
                cur_root = node_successor(cur_root);
            }
        }
        return self_t(_storage->_nodes.end(), _storage);
    }
    
    
    /**
     * @brief Return an iterator over all the subtrees for which `test(subtree)` returns
     * true. The iterator dereferences to `self_t`.
     *
     * The iterator's sequence finishes on (and excludes) this node's `subtree_end()` node.
     * `QueryIterator`s and `Subtree` nodes can be compared directly.
     *
     * `visit_all_nodes()` and `visit_all_leaf_nodes()` are convenience static member
     * functions of this class which can be passed to `test`.
     *
     * @param bound A callable object which accepts a `self_t` and returns `true` if 
     * that Subtree might contain objects which pass `test()`. Children of Subtrees
     * which do not pass `bound()` will be skipped.
     * @param test A callable object which accepts a `self_t` and returns `true` if
     * that Subtree should be visited by the iterator. If omitted, all nodes which 
     * pass `bound()` will be visited.
     */
    template <typename BoundingFn, typename TestFn>
    inline QueryIterator<self_t, BoundingFn, TestFn> query(BoundingFn bound, TestFn test=visit_all_nodes) const {
        return QueryIterator<self_t, BoundingFn, TestFn>(*this, bound, test);
    }
    
protected:
        
    // compute the position of the node that follows `node`'s subtree in 
    // "sibling-first" order. the successor to `node` occurs just before 
    // the next subtree to the right. an `O(k * log(n))` operation, on `n` 
    // the size of the tree, and `k` the arity. this is different from 
    // `++node` because `++node` may descend into the current node's 
    // subtree; we want to exactly hop over it.
    NodeRef node_successor(const NodeRef& node) const {
        NodeRef ancestor  = node->parent;
        NodeRef child     = node;
        // if there is no subtree to the right,
        // this is the conceptual successor:
        NodeRef successor = _storage->_nodes.end();
        
        while (ancestor != _storage->_nodes.end()) {
            NodeRef off_end_child = ancestor->child_last; 
            if (ancestor->n_children > 0) ++off_end_child;
            // look for a populated node to the right of us
            for (NodeRef c = child; c != off_end_child; ++c) {
                if (c->n_children > 0) {
                    successor = c->child_first;
                    break;
                }
            }
            if (successor != _storage->_nodes.end()) {
                // found the successor
                break;
            } else {
                // ascend the tree
                child    = ancestor;
                ancestor = ancestor->parent;
            }
        }
        return successor;
    }
    
    
    // compute the first item that occurs after the last item in the subtree
    // rooted at `node`, including if that subtree is empty.
    ItemRef item_successor(NodeRef node) const {
        for (NodeRef parent = node->parent; 
                parent != _storage->_nodes.end(); 
                node = parent, parent = node->parent) {
            NodeRef end_c = parent->child_last;
            if (parent->n_children != 0) ++end_c;
            for (NodeRef c = node; c != end_c; ++c) {
                if (c->n_items > 0) {
                    return c->items_first;
                }
            }
        }
        return this->_storage->_items.end();
    }
    
}; // class SubtreeBase



/**
 * @brief An const iterator to a subtree.
 */
template <typename NodeItem, typename LeafItem>
class ConstSubtree : public SubtreeBase<NodeItem, LeafItem, true> {
    
protected:
    
    typedef SubtreeBase<NodeItem, LeafItem, true> base_t;
    friend class SubtreeBase<NodeItem, LeafItem, true>;
    friend class Tree<NodeItem, LeafItem>;
    
    // allow myself to construct myself.
    // curiously recurring pattern is weird :[
    ConstSubtree(const base_t& other):base_t(other) {}
    
public:
    
    using typename base_t::iterator;
    using typename base_t::const_iterator;
    using typename base_t::item_iterator;
    using typename base_t::const_item_iterator;
    
    /// Construct a duplicate iterator to the same node of the same tree.
    ConstSubtree(const ConstSubtree<NodeItem, LeafItem>& other) = default;
    
    /// Construct a duplicate iterator to the same node of the same tree.
    ConstSubtree(const Subtree<NodeItem, LeafItem>& other):
        base_t(other._root, other._storage) {}
    
};



/**
 * @brief A non-const iterator to a subtree.
 */
template <typename NodeItem, typename LeafItem>
class Subtree : public SubtreeBase<NodeItem, LeafItem, false> {
    
protected:
    
    typedef SubtreeBase<NodeItem, LeafItem, false> base_t;
    using typename base_t::NodeRef;
    using typename base_t::ItemRef;
    using typename base_t::StorageRef;
    
    friend class SubtreeBase<NodeItem, LeafItem, false>;
    friend class Tree<NodeItem, LeafItem>;
    
    // allow myself to construct myself.
    Subtree(const base_t& other):base_t(other) {}
        
    Subtree(const NodeRef& root, StorageRef storage):
        base_t(root, storage) {}
        
public:
    
    typedef Subtree<NodeItem, LeafItem> self_t;
    using typename base_t::iterator;
    using typename base_t::const_iterator;
    using typename base_t::item_iterator;
    
    
    /// Construct a duplicate iterator to the same node of the same tree.
    Subtree(const Subtree<NodeItem, LeafItem>& other) = default;
    
    /// @name Mutation functions
    /// @{
    
    
    /**
     * @brief Insert a new tree node under this one, to the left of 
     * `insert_before`.
     *
     * If this node was previously empty, then its new first child will
     * contain all of its `LeafItem`s (this is to ensure those leaf
     * items all have a parent which is a leaf node). Otherwise, the new node
     * will be empty of any `LeafItem`s.
     *
     * If `insert_before` is the root node, or if `insert_before` is not 
     * a direct child (or end-child) of this node, then the tree is unaffected 
     * and `end()` is returned. Otherwise, the new node is returned.
     *
     * Potentially invalidates the `end()` of other subtrees, or the `begin()`
     * of this subtree if `insert_before` is `begin()`.
     *
     * @param insert_before A subtree pointing to a child or
     * end-child of this node.
     * @param args Construction arguments for the new `NodeItem`.
     * @return The newly created node, or `end()` if the node could not be 
     * created.
     */
    template <typename ... Args>
    self_t insert_child_node(
            const self_t& insert_before,
            Args&& ... args) const {
        NodeRef p = this->_root;
        NodeRef n = insert_before._root;
        bool parent_childless = p->n_children == 0;
        
        // check invalid cases
        if (n == this->_storage->_nodes.begin())       return this->end();
        if (n->parent != p and n != this->end()._root) return this->end();
        
        // set the insert point for empty nodes
        if (parent_childless) {
            n = this->node_successor(p);
        }
        
        // remember what our off-end node was
        NodeRef prev_end_node = this->end()._root;
        
        // create the new node
        NodeRef new_node = this->_storage->_nodes.insert(
            n, typename base_t::Node(this->_storage, p, args...));
        
        // take ownership of the parent's items, if we are its first child
        // (every item must belong to a leaf node, and the parent
        // has just become an internal node).
        if (parent_childless) {
            new_node->items_first = p->items_first;
            new_node->items_last  = p->items_last;
            new_node->n_items     = p->n_items;
        }
        
        // adjust child endpoints if necessary
        if (insert_before._root == prev_end_node) {
            p->child_last = new_node;
        }
        if (insert_before._root == p->child_first) {
            p->child_first = new_node;
        }
        p->n_children++;
        
        return self_t(new_node, this->_storage);
    }
    
    
    /**
     * @brief Convenience function to insert a new child node at the 
     * end of this node's children.
     *
     * @param obj NodeItem to insert into the new node.
     */
    inline self_t insert_child_node(const NodeItem& obj) const {
        return insert_child_node(this->end(), obj);
    }
    
    
    /**
     * @brief Insert new tree nodes to the left of `insert_before`, 
     * under this node. Return the first newly created node.
     *
     * If this node was previously empty, then its new first child will
     * contain all its `LeafItem`s (this is to ensure those leaf
     * items all have a parent which is a leaf node). All other new nodes
     * will be empty of any `LeafItem`s.
     *
     * If `insert_before` is the root node, or if `insert_before` is not 
     * a child (or end-child) of this node, then the tree is unaffected and 
     * `end()` is returned. Likewise, if no new nodes were inserted, 
     * `end()` is returned. Otherwise, the first new node is returned.
     * 
     * Potentially invalidates the `end()` of other subtrees, or the `begin()`
     * of this subtree if `insert_before` is `begin()`.
     *
     * @param insert_before A Subtree pointing to the sibling just 
     * after the last new node to be inserted.
     * @param i_begin A forward iterator to the first NodeItem to be inserted.
     * @param i_end A forward iterator just beyond the last NodeItem to be inserted.
     * @return The first newly created node, or `end()` if no new nodes were created.
     */
    template <typename NodeItemIterator>
    self_t insert_child_nodes(
            const self_t& insert_before,
                  NodeItemIterator  i_begin,
            const NodeItemIterator& i_end) const {
        NodeRef p = this->_root;
        NodeRef n = insert_before._root;
        bool parent_empty = p->n_children == 0;
        
        // check invalid cases
        if (n == this->_storage->_nodes.begin())       return this->end();
        if (n->parent != p and this->end()._root != n) return this->end();
        if (i_begin == i_end)                          return this->end();
        
        // set the insert point for empty nodes
        if (parent_empty) {
            n = this->node_successor(p);
        }
        
        // create the first node, and keep a reference to it
        // todo: this->_storage->_nodes.reserve(n_items);
        NodeRef prev_end_node  = this->end()._root;
        NodeRef first_new_node = this->_storage->_nodes.insert(
            n, base_t::Node(this->_storage, p, *(i_begin++)));
        NodeRef last_new_node  = first_new_node;
        
        // transfer ownership of the parent's items to the first child.
        // (every item must belong to a leaf node, and the parent
        // has just become an internal node).
        if (parent_empty) {
            first_new_node->items_first = p->items_first;
            first_new_node->items_last  = p->items_last;
            first_new_node->n_items     = p->n_items;
        }
        
        // add the rest of the nodes
        index_t ct = 1;
        for (NodeItemIterator i = i_begin; i != i_end; ++i, ++ct) {
            last_new_node = this->_storage->_nodes.insert(
                n, base_t::Node(this->_storage, p, *i));
        }
        
        // adjust child endpoints if necessary
        if (insert_before._root == prev_end_node) {
            p->child_last = last_new_node;
        }
        if (insert_before._root == p->child_first) {
            p->child_first = first_new_node;
        }
        p->n_children += ct;
        return self_t(first_new_node, this->_storage);
    }
    
    
    /** 
     * @brief Convenience function to insert multiple new child nodes at the
     * end of this node's children.
     */
    template <typename NodeItemIterator>
    inline self_t insert_child_nodes(
            const NodeItemIterator& i_begin,
            const NodeItemIterator& i_end) const {
        return insert_child_nodes(this->end(), i_begin, i_end);
    }
    
    
    /**
     * @brief Insert an item under this node.
     *
     * The new item will be inserted before the item at `insert_before`,
     * which must belong to this node, otherwise the behavior 
     * is undefined. It is permissible for `insert_before` to be this node's 
     * `items_begin()` or `items_end()`.
     * 
     * This must be a leaf node; i.e. one with no child nodes, so that 
     * placement of the new objects is not abiguous. If this is not a 
     * leaf node, the tree is unchanged.
     *
     * Potentially invalidates the `items_end()` of other subtrees, 
     * or the `items_begin()` of this subtree if `insert_before` is 
     * `items_begin()`.
     *
     * @param insert_before Item belonging to `parent` which the new items
     * are to be inserted before.
     * @param args New `LeafItem` to be inserted, or constructor arguments 
     * for the new `LeafItem`.
     *
     * @return An iterator to the first newly placed item if any were placed;
     * `items_end()` otherwise.
     */
    template <typename ... Args>
    item_iterator insert_item(
            const item_iterator& insert_before,
            Args&& ... args) const {
        // todo: it would be great if we could efficiently protect against 
        //   the user providing an insert_pt that doesn't belong to the parent
        
        // must insert to a leaf node
        if (this->node_count() > 0) {
            return this->items_end();
        }
        
        // find insert point
        NodeRef n = this->_root;
        ItemRef insert_pt = insert_before;
        if (n->n_items == 0) {
            insert_pt = this->item_successor(n);
        }
        
        // insert the item
        ItemRef new_item = this->_storage->_items.emplace(insert_pt, args...);
        
        // figure out the new begin/end items
        ItemRef new_begin, new_end;
        if (insert_before == n->items_first) {
            new_begin = new_item;
        } else {
            new_begin = n->items_first;
        }
        if (std::prev(new_item) == n->items_last or n->n_items == 0) {
            new_end = new_item;
        } else {
            new_end = n->items_last;
        }
        
        // udpate the boundaries of this node, and its ancestors
        this->_storage->update_boundary_items(
            this->_root,
            new_begin,
            new_end,
            1);
        
        return new_item;
    }
    
    
    /// Convenience function to insert a new item at the end.
    inline item_iterator insert_item(const LeafItem& item) const {
        return insert_item(this->items_end(), item);
    }
    
    
    /**
     * @brief Insert multiple leaf items into this node.
     *
     * The new items will be inserted before the item at `insert_before`, in 
     * the same order. `insert_before` must belong to the given node, 
     * otherwise the behavior is undefined. It is permissible for 
     * `insert_before` to be the node's `items_begin()` or `items_end()`.
     * 
     * This must be a leaf node; i.e. one with no child nodes, so that 
     * placement of the new objects is not abiguous. If this is not a 
     * leaf node, the tree is unchanged.
     *
     * Invalidates all `end()` `item_iterator`s.
     *
     * @param parent Leaf node under which the new items will be inserted.
     * @param insert_before Item belonging to `parent` which the new items
     * are to be inserted before.
     * @param first_item Forward iterator to first `LeafItem` to be inserted.
     * @param off_end_item Forward iterator just beyond the last `LeafItem` 
     * to be inserted.
     * @param new_item_count Optional return pointer to receive the count 
     * of newly placed objects.
     *
     * @return An iterator to the first newly placed item if any were placed;
     * `items_end()` otherwise.
     */
    template <typename LeafItemIterator>
    item_iterator insert_items(
            const item_iterator& insert_before,
            const LeafItemIterator  begin_item,
            const LeafItemIterator& end_item,
            index_t* new_item_count=nullptr) const {
        // todo: it would be great if we could protect against the user
        //   providing an insert_pt that doesn't belong to the parent
        
        // must insert to a leaf node
        if (this->node_count() > 0) {
            return this->items_end();
        }
        if (begin_item == end_item) {
            // nothing to be done
            return this->items_end();
        }
        
        // retrieve parent and insert point
        NodeRef n = this->_root;
        ItemRef insert_pt = insert_before;
        if (n.n_items == 0) {
            insert_pt = item_successor(n);
        }
        
        // insert the first item
        ItemRef first_new_item = this->_storage->_items.insert(insert_pt, *begin_item);
        ItemRef last_new_item  = first_new_item;
        ++begin_item;
        
        // insert the remaining items
        index_t item_ct = 1;
        for (LeafItemIterator i = begin_item; i != end_item; ++i) {
            last_new_item = this->_storage->_items.insert(insert_pt, *i);
            ++item_ct;
        }
        
        // figure out the new begin/end items
        ItemRef new_begin, new_end;
        if (insert_before == n->items_first) {
            new_begin = first_new_item;
        } else {
            new_begin = n->items_first;
        }
        if (std::prev(first_new_item) == n->items_last or n->n_items == 0) {
            new_end = last_new_item;
        } else {
            new_end = n->items_last;
        }
        
        // udpate the boundaries of this node, and its ancestors
        this->_storage->update_boundary_items(
            this->_root,
            new_begin,
            new_end,
            item_ct);
        
        if (new_item_count) *new_item_count = item_ct;
        return first_new_item;
    }
    
    
    /**
     * @brief Convenience function to insert multiple new items at
     * the end of this node's item list.
     */
    template <typename LeafItemIterator>
    inline item_iterator insert_items(
            const LeafItemIterator& begin_item,
            const LeafItemIterator& end_item) const {
        return insert_items(this->items_end(), begin_item, end_item);
    }
    
    
    /**
     * @brief Make a copy of `other_subtree` and insert it as a child of 
     * this node immediately before `insert_before`.
     *
     * If `insert_before` is not a child or end-child of this node, the
     * tree is unchanged and `end()` is returned.
     *
     * Invalidates all `end()` subtrees.
     *
     * @param other_subtree The subtree to copy and adopt as a child.
     * @param insert_before The child of this node before which to insert the new subtree.
     * @return An iterator pointing at the root of the new adopted subtree.
     */
    self_t adopt(const const_iterator& other_subtree, const self_t& insert_before) const {
        self_t new_node = insert_child_node(insert_before, *other_subtree);
        if (new_node == this->end()) return new_node;
        
        NodeRef n = new_node._root;
        
        // insert all the items
        // (also updates the boundary items; etc.)
        // new node is now a valid leaf with a buncha items:
        ItemRef first_new_item = n.insert_items(
            other_subtree.items_begin(),
            other_subtree.items_end());
        
        // todo: _storage->_nodes.reserve(...)
        
        // insert all the descendent nodes
        NodeRef subtree_end    = this->node_successor(n);
        NodeRef first_new_node = this->_storage->_nodes.insert(
            subtree_end, 
            other_subtree.subtree_begin()._root, 
            other_subtree.subtree_end()._root);
        
        // update the new root's child list
        n->n_children = other_subtree.node_count();
        if (first_new_node != subtree_end) {
            n->child_first = first_new_node;
            n->child_last  = subtree_end; --(n->child_last);
        }
        
        // all the pointers/iterators inside the subtree have to point
        // to *our* storage; not the original storage
        new_node.recalculate_references();
        
        return new_node;
    }
    
    
    /**
     * @brief Split this node into two sibling nodes.
     * 
     * A new node will be created just before this one, under the same parent, 
     * and this node's items will be split between them according to the 
     * result of `compare(item, pivot)`: If the comparison is less than zero, 
     * the item will be moved to the left (new) node. Otherwise, it will
     * remain in the right (existing) node.
     *
     * The split node must not be the root and must have no children.
     * If either is the case, the tree will not be changed and the off-end
     * node will be returned.
     *
     * All iterators to the items under this node are invalidated
     * by this operation, as well as potentially any `end()` `item_iterator`s.
     *
     * @param compare A callable object `compare(a,b)` which accepts a 
     * `LeafItem` as its left argument and a `P` as its right 
     * argument, and returns a signed number (negative for `a < b`, 
     * positive for `a > b`, zero for equality).
     * @param args Constructor arguments for the `NodeItem` object 
     * to be assigned to the newly-created (low) node.
     *
     * @return The newly created sibling, or the off-end node if this 
     *         is the root node or has child nodes.
     */
    template <typename Func, typename P, typename ... Args>
    self_t split(Func compare, P pivot, Args&& ... args) const {
        if (this->_root == this->_storage->_nodes.begin() or this->node_count() > 0) {
            return this->end();
        }
        
        self_t new_node = this->parent().insert_child_node(*this, args...);
        NodeRef hi = this->_root;
        NodeRef lo = new_node._root;
        
        // todo: if we had splice operations on list<>, we could avoid
        //       invalidating the item iterators.
        
        // traverse `node` looking for items to move to the new `lo` node.
        // because `lo` is adjacent to `hi`, we can add items to `lo` by moving the 
        // boundary between `lo` and `hi` and swapping the contents of the iterators.
        // we traverse backwards so that our swapping does not unsort our items.
        
        // corner case: hi is / becomes empty
        for (ItemRef i = hi->items_last; i != std::prev(hi->items_first); ) {
            if (compare(*i, pivot) < 0) {
                // move the item to the lower node
                
                // grab one item from the beginning of `hi`.
                // swap `*i` and the stolen item, so that the 
                // tested item lies in the low node where it belongs, 
                // and the stolen item lies under `i`. we'll check it
                // on the next iteration.
                
                if (lo->n_items == 0) {
                    lo->items_first = lo->items_last = hi->items_first;
                } else {
                    lo->items_last = hi->items_first;
                }
                
                ++(lo->n_items);
                --(hi->n_items);
                ++(hi->items_first);
                std::swap(*i, *(lo->items_last));
                
                // we do not decrement i, because we've just moved
                // an untested item underneath it. we should test
                // it first before moving on to the next one.
            } else {
                // `*i` is in the right place.
                --i;
            }
        }
        
        // if the `hi` node was emptied, fix up its item boundaries.
        if (hi->n_items == 0) {
            hi->items_first = hi->items_last = this->_storage->_items.end();
        }
        
        // we didn't delete/add any LeafItems, so we don't need to
        // update the boundary items.
        
        return new_node;
    }
    
    
    /**
     * @brief Remove an item from this subtree.
     *
     * If `item` is not in the subtree, or if this subtree is not a leaf
     * node, the tree will be unchanged.
     *
     * Invalidates the iterator to this item, and potentially any `end()`
     * `item_iterator`s.
     *
     * @param item The item to delete; a direct child of this node.
     * @param success Optional return value pointer to be filled with `true` 
     * if the item was deleted; `false` if the tree is unchanged.
     * @return An iterator pointing to the position of the item just beyond 
     * the one that was deleted. This will be `items_end()` if the 
     * deleted item was the last one in its parent node.
     */
    item_iterator erase(const item_iterator& item, bool* success=nullptr) const {
        item_iterator ret = this->items_end();
        bool ok = false;
        NodeRef n = this->_root;
        
        if (n->n_children == 0) {
            // verify that the item actually belongs to us.
            // (otherwise hell could break loose).
            // this is a linear check; if we *knew* we were the parent, the
            // entire operation would be log(n). laaaaame.
            bool found = false;
            for (item_iterator i = n->items_begin(); 
                    i != n->items_end(); 
                    ++i) {
                if (i == item) {
                    found = true;
                    break;
                }
            }
            
            if (found) {
                ItemRef next_item = this->_storage->_items.erase(item);
                ItemRef prev_item = next_item; --prev_item;
                
                // figure out the new begin / end items
                ItemRef new_begin, new_end;
                if (item == n->items_first) {
                    new_begin = next_item;
                } else {
                    new_begin = n->items_first;
                }
                if (item == n->items_last) {
                    new_end = prev_item;
                } else {
                    new_end = prev_item;
                }
                
                this->_storage->update_boundary_items(
                    n, 
                    new_begin,
                    new_end,
                    -1);
                
                ok  = true;
                ret = next_item;
            }
        }
        
        if (success) *success = ok;
        return ret;
    }
    
    
    /**
     * @brief Empty this node of all child nodes and descendent items.
     */
    void clear() const {
        flatten();
        
        NodeRef n = this->_root;
        
        if (n->n_items > 0) {
            // remember deleted references so they can be removed from ancestors
            ItemRef del_first_item = n->items_first;
            ItemRef del_last_item  = n->items_last;
            ItemRef endpt = n->items_last; ++endpt;
            
            // delete the items
            ItemRef next_item;
            index_t n_deleted;
            for (next_item = n->items_first; next_item != endpt; ++next_item) {
                this->_storage->_items.erase(next_item);
                ++n_deleted;
            }
            ItemRef prev_item = next_item; --prev_item;
            
            // update item boundaries / counts
            this->_storage->update_boundary_items(n, next_item, prev_item, n_deleted);
        }
    }
    
    
    /**
     * @brief Remove the subtree and its root at `child`, including all its leaf items.
     * 
     * If `child` is not a direct child of this node, the tree is unchanged.
     *
     * @param success An optional return value pointer; to be filled with `true`
     * if the child was deleted; `false` if the tree is unchanged.
     * @return An iterator pointing to the position of the node following the 
     * deleted one (which may be `end()`).
     */
    self_t erase(const self_t& child, bool* success=nullptr) const {
        // child must be our child
        if (child._root->parent != this->_root or 
            // and must not be global root
            child._root->parent == this->_storage->_nodes.end()) {
            if (success) *success = false;
            return this->end();
         }
        
        child.clear();
        
        NodeRef x = child._root;
        NodeRef next_node = this->_storage->_nodes.erase(x);
        NodeRef prev_node = next_node; --prev_node;
        
        // remove deleted child from its parent (us)
        this->_root->n_children--;
        if (this->_root->n_children == 0) {
            // `x` was the only child
            this->_root->child_first = 
                this->_root->child_last = 
                    this->_storage->_nodes.end();
        } else if (this->_root->child_first == x) {
            this->_root->child_first = next_node;
        } else if (this->_root->child_last == x) {
            this->_root->child_last = prev_node;
        }
        
        if (success) *success = true;
        return self_t(next_node, this->_storage);
    }
    
    
    /**
     * @brief Recursively remove all the children of this node, and take 
     * ownership of all items beneath them.
     */
    void flatten() const {
        if (this->node_count() == 0) return; // nothing to do.
        
        // because of sibling-first order, 
        // this is exactly all the descendents of `_root`:
        this->_storage->_nodes.erase(
            this->_root->child_first, 
            this->node_successor(this->_root));
        
        this->_root->n_children = 0;
        this->_root->child_first = 
            this->_root->child_last = 
                this->_storage->_nodes.end();
    }
    
    
    /// @}   // Mutation functions group
    
    
protected:
    
    // use all the child counts to recompute the first/last 
    // references in internal nodes. useful when taking ownership
    // of another tree's internal node list.
    // preconditions: 
    //   - first item AND child pointers of this node are correct
    //   - the descendents of this subtree have correct child and item counts
    //   - the nodes belonging to this subtree exist and are correctly ordered
    //     and placed in the tree
    //   - item boundaries of all ancestors are already correct
    //   - item counts of all ancestors are already correct
    // This should be an effective no-op if the subtree is already owned
    // (though it may still take O(n) to execute).
    void recalculate_references() const {
        NodeRef begin_child = this->_root->child_first;
        // next unassigned leaf item:
        ItemRef begin_items = this->_root->items_first;
        // last / most recent assigned leaf item in the list: 
        ItemRef last_item   = begin_items; --last_item;
        NodeRef end_child   = this->node_successor(this->_root);
        
        // in essence we "build" the tree by simply taking ownership of 
        // the nodes/items that already exist. we walk over the nodes in 
        // tree order, which means that child nodes are taken in order too.
        for (NodeRef node = this->_root; node != end_child; ) {
            node->child_first = node->child_last = 
                (node->n_children > 0) ? begin_child : this->_storage->_nodes.end();
            node->items_first = node->items_last = 
                (node->n_items    > 0) ? begin_items : this->_storage->_items.end();
            
            if (node->n_children > 0) {
                node->child_first->parent = node;
                
                // acquire direct children
                for (index_t c = 0; c < node->n_children - 1; ++c) {
                    ++(node->child_last);
                    node->child_last->parent = node;
                }
                
                // move the node "free space"
                begin_child = node->child_last; ++begin_child;
                
                // populate first child's subtree
                node = node->child_first;
            } else {
                if (node->n_items > 0) {
                    for (index_t i = 0; i < node->n_items - 1; ++i) {
                        ++(node->items_last);
                    }
                    // move the item "free space"
                    begin_items = last_item = node->items_last; ++begin_items;
                }
                    
                // fixup the items_end along the right edge of this subtree.
                NodeRef descendent = node;
                for (NodeRef ancestor = node->parent;
                        ancestor != this->_root;
                        descendent = ancestor, ancestor = ancestor->parent) {
                    if (ancestor->n_items > 0) ancestor->items_last = last_item;
                }
                
                // done populating this tree. populate the next sibling tree.
                node = this->node_successor(node);
            }
        }
    }
    
    
}; // ConstSubtree



/// @} // addtogroup storage


} // namespace geom
