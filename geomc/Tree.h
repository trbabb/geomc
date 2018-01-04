#include <geomc/geomc_defs.h>
#include <geomc/Templates.h>
#include <list>
#include <type_traits>

// todo: permit key/value structure.
//       nodes represent an abstract kind of set membership. node.contains(key). 
//       (maybe call it boundary so people don't retardedly nest data structures).
// todo: after key/value structure, permit a distance metric; allow KNN.
// todo: consider factoring the above simply as free functions which implement algorithms.
// todo: specialize for when leaf items or node items are void
// todo: smarter/more optimized underlying data structure
// todo: "flat array" variant tree. abstract NodeRef and ItemRef to be ptrs.
// todo: implement find(item_test(LeafItem) -> bool, node_test(NodeItem)=lambda x:true)
//       a generic algorithm for walking the tree to search for an item.
//       node_test() returns `true` if the node might possibly contain the item.
//       (templating over the functions would allow functors and inlining)
// todo: might be good to support re-rooting a tree at a child
//       or transplanting another part of the tree
// todo: might be good to support subtree-adoption from another tree
// todo: it is not possible to have an empty tree.

// xxx: recalculate_references needs to be factored for a non-root node.


// todo: make a function adopt(Tree& t)?
//       - does cheap splicing if owned by same tree
//       - otherwise does recalc_references on the subtree only.
//       - could replace copy by construction?


namespace geom {

    
/** @addtogroup storage
 *  @{
 */


// forward decls
template <typename NodeItem, typename LeafItem, bool Const> SubtreeBase;
template <typename NodeItem, typename LeafItem> Subtree;
template <typename NodeItem, typename LeafItem> ConstSubtree;


/**
 * @brief A dynamic tree of arbitrary arity.
 *
 * Each node in the tree may store an associated `NodeItem` object. 
 * A node may have any number of child nodes.
 * 
 * Leaf nodes may own a list of `LeafItem` objects. `LeafItem`s are kept
 * such that all the `LeafItems` in any subtree are contiguous, and internal nodes
 * provide iterators to the first and last `LeafItem`s in their respective subtrees.
 *
 * All mutation and access to a tree is provided by the Subtree and 
 * ConstSubtree classes; thin iterators pointing into the Tree data structure.
 *
 * @tparam NodeItem Type of data to be kept with internal tree nodes.
 * @tparam LeafItem Type of data to be kept by leaf nodes.
 */
template <typename NodeItem, typename LeafItem>
class Tree {
    
    template <bool Const> friend class ConstSubtree<NodeItem, LeafItem, Const>;
    friend class Subtree<NodeItem, LeafItem>;
    
    struct Node;
    
    typedef typename std::list<NodeItem>::iterator NodeRef;
    typedef typename std::list<LeafItem>::iterator ItemRef;
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
    
    
    // internally, nodes are kept in "sibling-first" order:
    //   - the successor of a node is its next sibling, if there is one,
    //   - otherwise the successor of a node is the first child of 
    //     the first sibling with a child, if there is one,
    //   - otherwise the successor is the first child of the first 
    //     subtree above and to the right, if there is one,
    //   - otherwise the successor is the end node.
    // as such, all direct children of a node are contiguous, and come 
    // after (though not necessarily *immediately* after) their parent.
    // all the descendents of a leftward child come before all the 
    // descendents of the child to its right. Therefore, any given subtree
    // (excluding its root) is contiguous. It is this recursive pattern:
    
    // [node] ... [c1 c2 c3 c4 ...]
    //    [children of c1] [descendents of c1...]
    //    [children of c2] [descendents of c2...]
    //    ...
    
    std::list<NodeItem> _nodes;
    std::list<LeafItem> _items;
    
    
    /// Construct an empty Tree.
    Tree() {
        _nodes.push_back(Node(this, _nodes.end()));
    }
    
    
    /// Default move constructor.
    Tree(const Tree<NodeItem, LeafItem>&& other) = default;
    
    
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
        NodeRef node = _nodes.front();
        
        // assign items to root
        node->items_first = _items.begin();
        node->items_last  = _items.end(); --(node->items_last);
        node->n_items     = m_nodes.size();
    }
    
    
    /// Make a copy of the given tree.
    Tree(const Tree<NodeItem, LeafItem>& other):
            m_nodes(other._nodes),
            m_items(other._items) {
        root().recalculate_references(); 
    }
    
    
    /// Construct a new tree from the given subtree.
    Tree(const ConstSubtree<NodeItem, LeafItem>& other):
            m_nodes(other.subtree_begin(), other.subtree_end()),
            m_items(other.items_begin(), other.items_end()) {
        // make a root (we only added the children)
        m_nodes.push_front(other._root);
        root().recalculate_references();
    }
    
    
    /// Get a subtree iterator to the root of this tree.
    Subtree<NodeItem, LeafItem> root() {
        return Subtree<NodeItem, LeafItem>(_nodes.front(), this);
    }
    
    
    /// Get a const subtree iterator to the root of this tree.
    ConstSubtree<NodeItem, LeafItem> root() const {
        return ConstSubtree<NodeItem, LeafItem>(_nodes.front(), this);
    }
    
        
    /// Tree equality
    bool operator==(const Tree<NodeItem, LeafItem>& other) const {
        if (this == &other) return true;
        if (item_count() != other.item_count()) return false;
        if (m_items      != other.m_items)      return false;
        
        // because the references will point to different memory,
        // we can't compare raw data. we have to examine the structure of the tree.
        NodeRef a_n =       m_nodes.begin();
        NodeRef b_n = other.m_nodes.begin();
        for (; a_n != m_nodes.end() and b_n != other.m_nodes.end(); ++a_n, ++b_n) {
            if (a_n->n_items    != b_n->n_items)    return false;
            if (a_n->n_children != b_n->n_children) return false;
            if (a_n->data       != b_n->data)       return false;
        }
        
        // must have same total number of nodes; not just a matching prefix subtree.
        // we don't check item_count() beforehand because it may not be stored, 
        // and will thus doubly incur the linear cost we have pay above.
        if (not (a_n == m_nodes.end() and b_n == other.m_nodes.end())) return false;
        
        return true;
    }
    
    
    /// Tree inequality
    bool operator!=(const Tree<NodeItem, LeafItem>& other) const {
        return not (*this == other);
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
    
    typedef Tree<NodeItem, LeafItem> Storage;
    typedef typename Storage::NodeRef       NodeRef;
    typedef typename Storage::ItemRef       ItemRef;
    typedef typename Storage::ItemRef       ConstItemRef;
    typedef typename Storage::Node          Node;
    typedef std::conditional<
            Const,
            const Storage*,
            Storage*>::type                 StorageRef;
    
    
    StorageRef _storage;
    NodeRef    _root;
    
    
    // construct an iterator to a particular node
    SubtreeBase(NodeRef& root, StorageRef storage):
        _storage(storage),
        _root(root) {}
    
public:
    
    /// Self type. A Subtree if this is a const iterator; a ConstSubtree otherwise.
    typedef typename std::conditional<
            Const,
            ConstSubtree<NodeItem, LeafItem, true>,
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
    typedef self_t          iterator;
    /// Const iterator over subtrees
    typedef ConstSubtree<NodeItem, LeafItem>  const_iterator;
    /// The tree's `NodeItem` type
    typedef NodeItem        value_type;
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
    
    
    /************************************
     * Methods                          *
     ************************************/

public:
    
    
    /// Obtain the tree into which this iterator points.
    inline typename std::conditional<Const, const tree_t&, tree_t&>::type
    tree() const {
        return _storage;
    }
    
    
    /// `+i`: Become first child
    inline self_t& operator+() {
        _root = _root->child_first;
        return *this;
    }
    
    /// `-i`: Become parent
    inline self_t& operator-() {
        _root = _root->parent;
        return *this;
    }
    
    /// `++i`: Become next sibling
    inline self_t& operator++() {
        _root = ++_root;
        return *this;
    }
    
    /// `--i`: Become previous sibling
    inline self_t& operator--() {
        _root = --_root;
        return *this;
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
    inline self_t subtree_begin() {
        return _root->child_first;
    }
    
    
    /**
     * @brief Return the last (off-end) subtree in the sequence of all subtrees
     * below this node.
     */
    inline self_t subtree_end() {
        return node_successor(_root);
    }
    
    
    /**
     * @brief Exhaustively search for the direct parent node of item `i` in 
     * this subtree.
     *
     * If `i` is not in the subtree under `node`, then return `end()`.
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
    
protected:   
        
    // compute the position of the node that follows `node`'s subtree in "sibling-first" order.
    // the successor to `node` occurs just before the next subtree to the right.
    // an `O(k * log(n))` operation, on `n` the size of the tree, and `k` the arity.
    // this is different from `++node` because `++node` may descend into the 
    // current node's subtree; we want to exactly hop over it.
    NodeRef node_successor(const NodeRef& node) const {
        NodeRef ancestor  = node->parent;
        NodeRef child     = node;
        NodeRef successor = nodes_end()._root; // if there is no subtree to the right,
                                               // this is the conceptual successor.
        while (ancestor != _storage->_nodes.end()) {
            NodeRef off_end_child = ancestor->child_last; 
            if (ancestor->n_children > 0) ++off_end_child;
            // look for a populated node to the right of us
            while (++child != off_end_child) {
                if (child->n_children > 0) {
                    successor = child->child_first;
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
            for (; node != end_c; ++node) {
                if (node->n_items > 0) {
                    return node->items_first;
                }
            }
        }
        return items_end();
    }
    
}; // class SubtreeBase



/**
 * @brief An const iterator to a subtree.
 */
template <typename NodeItem, typename LeafItem>
class ConstSubtree : public SubtreeBase<NodeItem, LeafItem, true> {
    
    /// Construct a duplicate iterator to the same node of the same tree.
    ConstSubtree(const ConstSubtree<NodeItem, LeafItem>& other) = default;
    
    
    /// Construct a duplicate iterator to the same node of the same tree.
    ConstSubtree(const Subtree<NodeItem, LeafItem& other>):
        SubtreeBase(other._root, other._storage) {}
    
};



/**
 * @brief A non-const iterator to a subtree.
 */
template <typename NodeItem, typename LeafItem>
class Subtree : public SubtreeBase<NodeItem, LeafItem, false> {
    
protected:
        
    Subtree(const NodeRef& root, StorageRef storage):
        Subtree(root, storage) {}
        
public:
    
    /// Construct a duplicate iterator to the same node of the same tree.
    Subtree(const Subtree<NodeItem, LeafItem>& other) = default;
    
    
    /**
     * @brief Insert an item under this node.
     *
     * The new item will be inserted before the item at `insert_before`.
     * `insert_before` must belong to this node, otherwise the behavior 
     * is undefined. It is permissible for `insert_before` to be this node's 
     * `items_begin()` or `items_end()`.
     * 
     * This must be a leaf node; i.e. one with no child nodes, so that 
     * placement of the new objects is not abiguous. If this is not a 
     * leaf node, the tree is unchanged.
     *
     * Invalidates all `end()` `node_iterator`s.
     *
     * @param insert_before Item belonging to `parent` which the new items
     * are to be inserted before.
     * @param obj New leaf item to be inserted.
     *
     * @return An iterator to the first newly placed item if any were placed;
     * `items_end()` otherwise.
     */
    item_iterator insert_item(
            const item_iterator& insert_before,
            const LeafItem& obj) {
        // todo: it would be great if we could protect against the user
        //   providing an insert_pt that doesn't belong to the parent
        
        // must insert to a leaf node
        if (node_count() > 0) {
            return items_end();
        }
        
        // find insert point
        NodeRef n = _root;
        ItemRef insert_pt = insert_before;
        if (n->n_items == 0) {
            insert_pt = item_successor(n);
        }
        
        // insert the item
        ItemRef new_item = _storage->_nodes.insert(insert_pt, obj);
        
        // ascend the tree, updating the item ranges of the ancestors
        ItemRef precursor_item = new_item; --precursor_item;
        while (n != _storage->_nodes.end()) {
            if (n->items_first == insert_pt or n->n_items == 0) {
                n->items_first = new_item;
            }
            if (n->items_last == precursor_item or n->n_items == 0) {
                n->items_last = new_item;
            }
            n->n_items += 1;
            n = n->parent;
        }
        
        return new_item;
    }
    
    
    /// Convenience function to insert a new item at the end.
    inline item_iterator insert_item(const LeafItem& item) {
        return insert_item(items_end(), item);
    }
    
    
    /**
     * @brief Insert multiple leaf items into this node.
     *
     * The new items will be inserted before the item at `insert_before`, in 
     * the same order. `insert_before` must belong to the given node, 
     * otherwise the behavior is undefined. It is permissible for 
     * `insert_before` to be the node's `items_begin()` or `items_end()`.
     * 
     * `node` must be a leaf node; i.e. one with no child nodes, so that 
     * placement of the new objects is not abiguous. If `node` is not a 
     * leaf node, the tree is unchanged.
     *
     * Invalidates all `end()` `node_iterator`s.
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
            index_t* new_item_count=nullptr) {
        // todo: it would be great if we could protect against the user
        //   providing an insert_pt that doesn't belong to the parent
        
        // must insert to a leaf node
        if (node_count() > 0) {
            return items_end();
        }
        if (begin_item == end_item) {
            // nothing to be done
            return items_end();
        }
        
        // retrieve parent and insert point
        NodeRef n = _root;
        ItemRef insert_pt = insert_before;
        if (n.n_items == 0) {
            insert_pt = item_successor(n);
        }
        
        // insert the first item
        ItemRef first_new_item = _storage->_items.insert(insert_pt, *begin_item);
        ItemRef last_new_item  = first_new_item;
        ++begin_item;
        
        // insert the remaining items
        index_t item_ct = 1;
        for (LeafItemIterator i = begin_item; i != end_item; ++i) {
            last_new_item = _storage->_items.insert(insert_pt, *i);
            ++item_ct;
        }
        
        // ascend the tree, updating the item ranges of the ancestors
        ItemRef precursor_item = first_new_item; --precursor_item;
        while (n != _storage->_nodes.end()) {
            if (n->items_first == insert_pt or n->n_items == 0) {
                n->items_first = first_new_item;
            }
            if (n->items_last == precursor_item or n->n_items == 0) {
                n->items_last = last_new_item;
            }
            n->n_items += item_ct;
            n = n->parent;
        }
        
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
            const LeafItemIterator& end_item) {
        return insert_items(items_end(), begin_item, end_item);
    }
    
    
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
     * Invalidates all `end()` `node_iterator`s.
     *
     * @param insert_before A `node_iterator` pointing to a child or
     * end-child of this node.
     * @param args Construction arguments for the new `NodeItem`.
     * @return The newly created node, or `end()` if the node could not be 
     * created.
     */
    template <typename ... Args>
    self_t insert_child_node(
            const self_t& insert_before,
            Args&& ... args) {
        NodeRef p = _root;
        NodeRef n = insert_before._root;
        bool parent_childless = p->n_children == 0;
        
        // check invalid cases
        if (n == _storage->_nodes.begin())       return self_t(_storage->_nodes.end(), _storage);
        if (n->parent != p and n != end()._root) return self_t(_storage->_nodes.end(), _storage);
        
        // set the insert point for empty nodes
        if (parent_childless) {
            n = node_successor(p);
        }
        
        // remember what our off-end node was
        NodeRef prev_end_node = end()._root;
        
        // create the new node
        NodeRef new_node = _storage->_nodes.insert(n, Node(_storage, p, args...));
        
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
        
        return self_t(new_node, _storage);
    }
    
    
    /**
     * @brief Convenience function to insert a new child node.
     */
    inline self_t insert_child_node(const NodeItem& obj) {
        return insert_child_node(end(), obj);
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
     * Invalidates all `end()` Subtrees.
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
            const NodeItemIterator& i_end) {
        NodeRef p = _root;
        NodeRef n = insert_before._root;
        bool parent_empty = p->n_children == 0;
        
        // check invalid cases
        if (n == _storage->_nodes.begin())       return self_t(_storage->_nodes.end(), _storage);
        if (n->parent != p and end()._root != n) return self_t(_storage->_nodes.end(), _storage);
        if (i_begin == i_end)                    return self_t(_storage->_nodes.end(), _storage);
        
        // set the insert point for empty nodes
        if (parent_empty) {
            n = node_successor(p);
        }
        
        // create the first node, and keep a reference to it
        // todo: _nodes.reserve(n_items);
        NodeRef prev_end_node  = end()._root;
        NodeRef first_new_node = _storage->_nodes.insert(n, Node(_storage, p, *(i_begin++)));
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
            last_new_node = _storage->_nodes.insert(n, Node(_storage, p, *i));
        }
        
        // adjust child endpoints if necessary
        if (insert_before._root == prev_end_node) {
            p->child_last = last_new_node;
        }
        if (insert_before._root == p->child_first) {
            p->child_first = first_new_node;
        }
        p->n_children += ct;
        return self_t(first_new_node, _storage);
    }
    
    
    /** 
     * @brief Convenience function to insert multiple new child nodes after
     * the last node under their parent.
     */
    template <typename NodeItemIterator>
    inline self_t insert_child_nodes(
            const NodeItemIterator& i_begin,
            const NodeItemIterator& i_end) {
        return insert_child_nodes(end(), i_begin, i_end);
    }
    
    
    /**
     * @brief Split `node` into two sibling nodes.
     * 
     * A new node will be created just before `node` under the same parent, 
     * and the items of `node` will be split between them according to the 
     * result of `compare(item, pivot)`: If the comparison is less than zero, 
     * the item will be moved to the left (new) node. Otherwise, it will
     * remain in the right (existing) node.
     *
     * The split node must not be the root and must have no children.
     * If either is the case, the tree will not be changed and the off-end
     * node will be returned.
     *
     * All iterators to the items under `node` are invalidated
     * by this operation, as well as all `end()` `item_iterator`s.
     *
     * @param compare A callable object `compare(a,b)` which accepts a 
     *        `LeafItem` as its left argument and a `P` as its right 
     *        argument, and returns a signed number (negative for `a < b`, 
     *        positive for `a > b`, zero for equality).
     * @param obj The internal node object to assign to the newly-created 
     *        (low) node.
     *
     * @return The newly created sibling, or the off-end node if `node` 
     *         is the root node or has child nodes.
     */
    template <typename Func, typename P>
    self_t split(Func compare, P pivot, NodeItem& obj) {
        if (_root == _storage->_nodes.begin() or node_count() > 0) {
            return self_t(_storage->_nodes.end(), _storage);
        }
        
        self_t new_node = insert_child_node(_root->parent, _root, obj);
        NodeRef hi = _root;
        NodeRef lo = new_node._root;
        
        // todo: if we had splice operations on list<>, we could avoid
        //       invalidating the item iterators.
        
        // traverse `node` looking for items to move to the new `lo` node.
        // because `lo` is adjacent to `hi`, we can add items to `lo` by moving the 
        // boundary between `lo` and `hi` and swapping the contents of the iterators.
        // we traverse backwards so that our swapping does not unsort our items.
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
                
                ++lo->n_items;
                --hi->n_items;
                ++hi->items_first;
                std::swap(*i, *(lo->items_last));
                
                // we do not decrement i, because we've just moved
                // an untested item underneath it. we should test
                // it first before moving on to the next one.
            } else {
                // `*i` is in the right place.
                --i;
            }
        }
        return new_node;
    }
    
    
    /**
     * @brief Remove an item from this subtree.
     *
     * If `item` is not in the subtree, the tree will be unchanged.
     *
     * Calling this method on a close ancestor of `item` reduces the 
     * cost of finding `item`'s parent in the tree before removal. 
     * In general, an `O(n * log(n))` search is performed to find the 
     * direct parent of `item`, where `n` is the size of this subtree.
     *
     * Invalidates the iterator to this item, as well as any `end()`
     * `item_iterator`s.
     *
     * @return An iterator pointing to the position of the item just beyond 
     * the one that was deleted. This will be `parent.items_end()` if the 
     * deleted item was the last one in its parent node.
     */
    item_iterator erase(const item_iterator& item) {
        NodeRef parent = find_parent(item)._root;
        
        if (parent != _storage->_nodes.end()) return _erase(item, parent);
        
        return _storage->_items.end();
        
        // todo: in key/value subclass, you can find `item` in log(n) time. 
        // be sure to override this fucker.
        
        // todo: no way to detect error condition of "item not found".
    }
    
    
    /**
     * @brief Empty this node of all child nodes and descendent items.
     */
    void clear() {
        flatten();
        
        NodeRef n = _root;
        
        if (n->n_items > 0) {
            // remember deleted references so they can be removed from ancestors
            ItemRef del_first_item = n->items_first;
            ItemRef del_last_item  = n->items_last;
            ItemRef endpt = n->items_last; ++endpt;
            
            // delete the items
            ItemRef next_item;
            index_t n_deleted;
            for (next_item = n->items_first; next_item != endpt; ++next_item) {
                _storage->_items.erase(next_item);
                ++n_deleted;
            }
            ItemRef prev_item = next_item; --prev_item;
            
            // update ancestors' item count and boundary items
            while (n != _storage->_nodes.end()) {
                n->n_items -= n_deleted;
                // fix up the boundary items
                if (n->n_items == 0) {
                    n->items_first = n->items_last = _storage->_items.end();
                } else {
                    if (n->items_first == del_first_item) {
                        n->items_first = next_item;
                    }
                    if (n->items_last == del_last_item) {
                        n->items_last = prev_item;
                    }
                }
                n = n->parent;
            }
        }
    }
    
    
    /**
     * @brief Remove the subtree and its root at `child`, including all its leaf items.
     * 
     * If `child` is not a direct child of this node, the tree is unchanged.
     *
     * @return An iterator pointing to the position of the node following the 
     * deleted one (which may be `parent.end()`).
     */
    self_t erase(self_t& child) {
        // child must be our child
        if (child._root->parent != _root) {
            return self_t(_storage->_nodes.end(), _storage);
        }
        // child must not be global root
        if (child._root->parent == _storage->_nodes.end()) {
            return self_t(_storage->_nodes.end(), _storage);
        }
        
        child.clear();
        
        NodeRef x = child._root;
        NodeRef next_node = _storage->_nodes.erase(x);
        NodeRef prev_node = next_node; --prev_node;
        
        // remove deleted child from its parent (us)
        _root->n_children--;
        if (_root->n_children == 0) {
            // `x` was the only child
            _root->child_first = n->child_last = _storage->_nodes.end();
        } else if (_root->child_first == x) {
            _root->child_first = next_node;
        } else if (_root->child_last == x) {
            _root->child_last = prev_node;
        }
        
        return self_t(next_node, _storage);
    }
    
    
    /**
     * @brief Recursively remove all the children of this node, and take 
     * ownership of all items beneath them.
     */
    void flatten() {
        if (node_count() == 0) return; // nothing to do.
        
        // because of sibling-first order, 
        // this is exactly all the descendents of `_root`:
        _storage->_nodes.erase(_root->child_first, node_successor(_root));
        
        _root->n_children = 0;
        _root->child_first = _root->child_last = _storage->_nodes.end();
    }
    
    
protected:
    
    // erase when `parent` is known to be the exact parent of `item`; not just an ancestor.
    item_iterator _erase(const ItemRef& item, NodeRef parent) {
        ItemRef next_item = _storage->_items.erase(item);
        ItemRef prev_item = next_item; --prev_item;
        
        // update the ancestors.
        while (parent != _storage->_nodes.end()) {
            parent->n_items--;
            if (parent->n_items == 0) {
                parent->items_first = parent->items_last = _storage->_items.end();
            } else {
                if (parent->child_first == item) {
                    parent->child_first = next_item;
                }
                if (parent->child_last == item) {
                    parent->child_last = prev_item;
                }
            }
            parent = parent->parent;
        }
        
        return next_item;
    }
    
    // todo: refactor this to work on an internal node
    
    // use all the child counts to recompute the first/last 
    // references in internal nodes. useful when taking ownership
    // of another tree's internal node list.
    void recalculate_references() {
        NodeRef nod = _storage->_nodes.begin();
        NodeRef begin_child = nod; ++begin_child;
        ItemRef begin_items = _storage->_items.begin(); // next unassigned leaf item
        ItemRef last_item   = _storage->_items.end();   // last / most recent assigned leaf item in the list
        
        // in essence we "build" the tree by simply taking ownership of 
        // the nodes/items that already exist. we walk over the nodes in 
        // tree order, which means that child nodes are taken in order too.
        // this spares us recursion (i.e. no possibility of stack-blowing) 
        // and some logic keeping track of descent/ascent.
        for (; nod != _storage->_nodes.end(); ++nod) {
            nod->child_first = nod->child_last = (nod->n_children > 0) ? begin_child : _storage->_nodes.end();
            nod->items_first = nod->items_last = (nod->n_items    > 0) ? begin_items : _storage->_items.end();
            
            if (nod->n_children > 0) {
                nod->child_first->parent = nod;
                
                // acquire direct children
                for (index_t c = 0; c < nod->n_children - 1; ++c) {
                    ++nod->child_last;
                    nod->child_last->parent = nod;
                }
                
                // move the node "free space"
                begin_child = nod->child_last; ++begin_child;
            } else {
                if (nod->n_items > 0) {
                    for (index_t i = 0; i < nod->n_items - 1; ++i) {
                        ++nod->items_last;
                    }
                    // move the item "free space"
                    begin_items = last_item = nod->items_last; ++begin_items;
                }
                    
                // fixup the items_end along the right edge of this subtree.
                NodeRef descendent = nod;
                for (NodeRef ancestor = nod->parent; 
                        ancestor != _storage->_nodes.end(); 
                        descendent = ancestor, ancestor = ancestor->parent) {
                    if (ancestor->child_last == descendent) {
                        if (ancestor->n_items > 0) ancestor->items_last = last_item;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    
    
}; // ConstSubtree




/// @} // addtogroup storage


} // namespace geom
