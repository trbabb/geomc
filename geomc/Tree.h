#include <geomc/Templates.h>

// todo: permit key/value structure.
//       nodes represent an abstract kind of set membership. node.contains(key). 
//       (maybe call it boundary so people don't retardedly nest data structures).
// todo: after key/value structure, permit a distance metric; allow KNN.
// todo: specialize for when leaf items or node items are void
// todo: smarter/more optimized underlying data structure
// todo: insert multiple items
// todo: "flat array" variant tree. abstract NodeRef and ItemRef to be ptrs.
// todo: figure out a way to put node_successor, aka subtree_last into node iterator?

namespace geom {
    
/** @addtogroup storage
 *  @{
 */

/**
 * @brief A dynamic tree of arbitrary arity.
 *
 * Each node in the tree may store an associated `NodeItem` object. A node may have any number of child nodes.
 * 
 * Leaf nodes may own a list of `LeafItem` objects. `LeafItem`s are kept
 * such that all the `LeafItems` in a subtree are contiguous, and internal nodes
 * can provide iterators to the first and last `LeafItem`s in their respective subtrees.
 *
 * @tparam NodeItem Type of data to be kept with internal tree nodes.
 * @tparam LeafItem Type of data to be kept by leaf nodes.
 */
template <typename NodeItem, typename LeafItem>
class Tree {
    
    // fwd decl
    struct Node;
    
    typedef typename std::list<Node>::iterator     NodeRef;
    typedef typename std::list<LeafItem>::iterator ItemRef;
    
    struct Node {
        
        Node() {}
        
        Node(Tree<NodeItem, LeafItem>* t,
             NodeRef parent,
             const NodeItem& obj): // xxx move semantics?
                parent(parent),
                child_first (t->m_nodes.end()),
                child_last  (t->m_nodes.end()), 
                items_first (t->m_items.end()),
                items_last  (t->m_items.end()),
                n_items     (0),
                n_children  (0),
                data        (obj) {}
        
        NodeRef  parent;
        NodeRef  child_first;
        NodeRef  child_last;
        ItemRef  items_first;
        ItemRef  items_last;
        
        index_t  n_items;
        index_t  n_children;
        
        NodeItem data;
    };
    
    
    // internally, nodes are kept in (sort of) breadth-first order.
    // all direct children of a node are contiguous, and come 
    // after (though not necessarily *immediately* after) their parent.
    // all the descendents of a leftward child come before all the 
    // descendents of the child to its right. Therefore, any given subtree
    // (excluding its root) is contiguous. It is this recursive pattern:
    
    // [node] [c1 c2 c3 c4 ...] [children of c1] [descendents of c1...] [children of c2] [descendents of c2...] ...
      
    std::list<Node>     nodes;
    
    // items are stored in depth-first order. 
    // all items in a subtree are contiguous.
    
    std::list<LeafItem> items;
    
public:
    
    template <bool Const>
    class ItemIterator {
        
        friend class Tree<NodeItem,LeafItem>;
        
        typedef ItemIterator<Const> self_t;
        typedef typename ConstType<LeafItem,Const>::reference_t reference; // reference to value type
        typedef typename ConstType<LeafItem,Const>::pointer_t   pointer;   // pointer to value type
        
        ItemRef item;
        NodeRef parent; // an ancestor of the item. may or may not be "close" depending on how the item was found.
        
        ItemIterator(const ItemRef& item, const NodeRef& parent):
            item(item),
            parent(parent) {}
        
        // can always promote to const:
        ItemIterator(const ItemIterator<false>& other):
            item(other.item),
            parent(other.parent) {}
        
    public:
        
        inline self_t& operator++() {
            ++item;
            return *this;
        }
        
        inline self_t& operator--() {
            --item;
            return *this;
        }
        
        inline self_t operator++(int) {
            self_t tmp = item;
            ++item;
            return tmp;
        }
        
        inline self_t operator--(int) {
            self_t tmp = item;
            --item;
            return tmp;
        }
        
        inline reference operator*() const {
            return *item;
        }
        
        inline pointer operator->() const {
            return &(*item);
        }
        
        template <bool C>
        inline bool operator==(const ItemIterator<C>& other) const {
            return item == other.item;
        }
        
        template <bool C>
        inline bool operator!=(const ItemIterator<C>& other) const {
            return item != other.item;
        }
        
    };
    
    
    /**
     * @brief An optionally-const iterator over the internal nodes of a tree.
     *
     * Dereferencing this iterator produces a `NodeItem` object.
     *
     * This class can be treated as a container of objects of its own type.
     *
     * @tparam Const Whether this iterator refers to a `const` item or not.
     */
    template <bool Const>
    class NodeIterator {
        
        friend class Tree<NodeItem,LeafItem>;
        
        NodeRef node;
        
        NodeIterator(NodeRef n) : node(n) {}
        
        // can always promote to const:
        NodeIterator(const NodeIterator<false>& i) : node(i.node) {}
        
    public:
        
        /// The tree's `NodeItem` type.
        typedef NodeItem                                         value_type;
        /// A (possibly const) reference to a `NodeItem`.
        typedef typename ConstType<NodeItem,Const>::reference_t  reference;       // reference to value type
        /// A (possibly const) pointer to a `NodeItem`.
        typedef typename ConstType<NodeItem,Const>::pointer_t    pointer;         // pointer to value type
        /// A const reference to a `NodeItem`.
        typedef const NodeItem&                                  const_reference; // const ref to value type
        /// Iterator over child nodes; same as self type.
        typedef NodeIterator<Const>                              iterator;        // self type (container concept)
        
        typedef NodeIterator<Const>                              self_t;          // self type
        
        /// Iterator over `LeafItem`s.
        typedef ItemIterator<const> item_iterator; // iterator for leaf items
        
        /// `+i`: Become first child
        inline self_t& operator+() {
            node = node->child_first;
            return *this;
        }
        
        /// `-i`: Become parent
        inline self_t& operator-() {
            node = node->parent;
            return *this;
        }
        
        /// `++i`: Become next sibling
        inline self_t& operator++() {
            node = ++node;
            return *this;
        }
        
        /// `--i`: Become previous sibling
        inline self_t& operator--() {
            node = --node;
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
            NodeIterator tmp = *this;
            --(*this);
            return tmp;
        }
        
        /// `*i`: Get node value
        inline reference operator*() const {
            return node->data;
        }
        
        /// `i->...`: Access node value member
        inline pointer operator->() const {
            return &node->data;
        }
        
        /// Number of direct child nodes
        inline index_t node_count() const {
            return node->n_children;
        }
        
        /// Number of leaf items in this subtree
        inline index_t item_count() const {
            return node->n_items;
        }
        
        /// Get first child
        inline self_t begin() const {
            return node->child_first;
        }
        
        /// Get last (off-end) child
        inline self_t end() const {
            NodeRef tmp = node->child_last;
            return ++tmp;
        }
        
        /// Get first object inside this subtree
        inline item_iterator items_begin() const {
            return item_iterator(node->items_first, node);
        }
        
        /**
         * @brief Get last (off-end) object in this subtree. It is invalid to increment or dereference this iterator.
         */
        inline item_iterator items_end() const {
            item_iterator tmp = node->items_last;
            return item_iterator(++tmp, node);
        }
    };
    
    
    /// Iterator over the internal tree nodes
    typedef NodeIterator<false> node_iterator;
    
    /// Const iterator over the internal tree nodes
    typedef NodeIterator<true>  const_node_iterator;
    
    /// Iterator over `LeafItem`s
    typedef ItemIterator<false> item_iterator;
    
    /// Const iterator over `LeafItem`s
    typedef ItemIterator<true> const_item_iterator;
    
    
    typedef Tree<NodeItem, LeafItem> self_t;
    
    
private:
    
    /************************************
     * Members                          *
     ************************************/
    
    // invariant after construction: root node exists.
    
    std::list<Node>     m_nodes;
    std::list<LeafItem> m_items;

public:
    
    
    /// Construct an empty Tree.
    Tree() {
        nodes.push_back(Node());
        Node&       n = m_nodes.front();
        n.parent      = m_nodes.end();
        n.child_first = m_nodes.end(); // no children
        n.child_last  = m_nodes.end();
        n.items_first = m_items.end(); // no items
        n.items_last  = m_items.end();
        n.n_items     = 0;
        n.n_children  = 0;
    }
    
    
    /**
     * @brief Construct a tree having a single root node, and copy the
     * items in `[begin, end)` to it.
     */
    template <typename LeafItemIterator>
    Tree(LeafItemIterator begin, LeafItemIterator end) {
        m_items.assign(begin, end);
        
        nodes.push_back(Node());
        Node&       n = m_nodes.front();
        n.parent      = m_nodes.end();
        n.child_first = m_nodes.end();
        n.child_last  = m_nodes.end();
        n.items_first = m_items.begin();
        n.items_last  = m_items.end(); --n.items_last;
        n.n_items     = m_nodes.size();
        n.n_children  = 0;
    }
    
    
    /// Make a copy of the given tree.
    Tree(const Tree<NodeItem, LeafItem>& other):
            m_nodes(other.m_nodes),
            m_items(other.m_items) {
        recalculate_references(); 
    }
    
    
    /// Construct a new tree from the subtree in `tree` rooted at `n`.
    Tree(const Tree<NodeItem, LeafItem>& tree, node_iterator n):
            m_nodes(n.begin().node, tree.node_successor(n.node)),
            m_items(n.items_begin().item, n.items_end().item) {
        m_nodes.push_front(n.node);
        recalculate_references();
    }
    
    
    /// Number of leaf items in the tree.
    inline index_t item_count() const {
        return m_nodes.begin()->n_items;
    }
    
    
    /// Number of nodes in the tree.
    inline index_t node_count() const {
        // xxx: may be slow for stupid implementations of std::list. :(
        // we can't use root->n_children, because that's only *direct* children.
        return m_nodes.size();
    }
    
    
    /**
     * @brief Return an iterator pointing at the root tree node. 
     * 
     * Use `+i` and `-i` to ascend and descend to the first-child and parent nodes
     * respectively. `++i` and `--i` navigate the tree in breadth-first order, and can
     * therefore be used to find the next and previous siblings.
     */
    inline node_iterator nodes_begin() {
        return m_nodes.begin();
    }
    
    
    /**
     * @brief Return an iterator pointing beyond the last node in the tree.
     * 
     * This node is conceptually the parent of the root.
     */
    inline node_iterator nodes_end() {
        return m_nodes.end();
    }
    
    
    /**
     * @brief Return a const iterator pointing at the root tree node. 
     * 
     * Use `+i` and `-i` to ascend and descend to the first-child and parent nodes
     * respectively. `++i` and `--i` navigate the tree in breadth-first order, and can
     * therefore be used to find the next and previous siblings.
     */
    inline const_node_iterator nodes_begin() const {
        return const_cast<self_t*>(this)->m_nodes.begin();
    }
    
    
    /**
     * @brief Return an iterator pointing beyond the last node in the tree.
     * 
     * This node is conceptually the parent of the root.
     */
    inline const_node_iterator nodes_end() const {
        return const_cast<self_t*>(this)->m_nodes.end();
    }
    
    
    /**
     * @brief Return an iterator pointing at the first leaf item in the tree.
     */
    inline item_iterator items_begin() {
        return item_iterator(m_items.begin(), m_nodes.begin());
    }
    
    
    /**
     * @brief Return an iterator just beyond the last leaf item in the tree.
     */
    inline item_iterator items_end() {
        return item_iterator(m_items.end(), m_nodes.begin());
    }
    
    
    /**
     * @brief Return an iterator pointing at the first (const) leaf item in the tree.
     */
    inline const_item_iterator items_begin() const {
        self_t* t = const_cast<self_t*>(this);
        return const_item_iterator(
            t->m_items.begin(), 
            t->m_nodes.begin());
    }
    
    
    /**
     * @brief Return an iterator just beyond the last (const) leaf item in the tree.
     */
    inline const_item_iterator items_end() const {
        self_t* t = const_cast<self_t*>(this);
        return const_item_iterator(
            t->m_items.end(), 
            t->m_nodes.begin());
    }
    
    
    /**
     * @brief Return the first node in a contiguous block of nodes in the subtree rooted at `root` 
     * (excluding `root` itself).
     */
    inline node_iterator subtree_begin(const node_iterator& root) {
        return root->child_first;
    }
    
    
    /**
     * @brief Return the last (off-end) node in a contiguous block of nodes in the subtree rooted at `root`.
     */
    inline node_iterator subtree_end(const node_iterator& root) {
        return node_successor(root.node);
    }
    
    
    /**
     * @brief Return the first node in a contiguous block of nodes in the subtree rooted at `root` 
     * (excluding `root` itself).
     */
    inline const_node_iterator subtree_begin(const node_iterator& root) const {
        return root->child_first;
    }
    
    
    /**
     * @brief Return the last (off-end) node in a contiguous block of nodes in the subtree rooted at `root`.
     */
    inline const_node_iterator subtree_end(const node_iterator& root) const {
        return node_successor(root.node);
    }
    
    
    /**
     * @brief Exhaustively search for the direct parent node of item `i` in the subtree rooted at `ancestor`.
     *
     * If `i` is not in the subtree under `node`, then return `end()`.
     */
    node_iterator find_parent(const node_iterator& ancestor, const item_iterator& i) const {
        NodeRef subtree_last = node_successor(ancestor.node);
        for (NodeRef n = ancestor.node; n != subtree_last; ++n) {
            if (n->n_children == 0) {
                ItemRef end_item = n->items_last; ++end_item;
                for (ItemRef j = n->items_first; j != end_item; ++j) {
                    if (j == i.item) return n;
                }
            }
        }
        return m_nodes.end();
    }
    
    
    /**
     * @brief Insert the given item into the tree under the given node.
     * 
     * `node` must be a leaf node; i.e. one with no child nodes, so that placement of `obj` is not ambiguous. 
     * If `node` is not a leaf node, the tree is unchanged.
     *
     * @return An iterator to the newly placed item, or `items_end()` if `node` is not a leaf node.
     */
    item_iterator insert(const node_iterator& parent, const LeafItem& obj) { // todo: move semantics?
        if (parent.nodes() > 0) return items.end();
        
        NodeRef node = parent.node;
        ItemRef prev_item;
        ItemRef insert_pos;
        
        if (parent.item_count() > 0) {
            // node already has items.
            prev_item = node->items_last;
            insert_pos = prev_item; ++insert_pos;
        } else {
            // adding the first item to an empty node.
            insert_pos = item_successor(node);
            prev_item = insert_pos; --prev_item;
        }
        
        ItemRef new_item = m_items.insert(insert_pos, obj);
        
        // walk upwards from the item's parent, updating the item count and boundary items.
        while (node != m_nodes.end()) {
            if (node->items_first == insert_pos or node->n_items == 0) {
                node->items_first = new_item;
            }
            if (node->items_last == prev_item or node->n_items == 0) {
                node->items_last = new_item;
            }
            node->n_items++;
            node = node->parent;
        }
        
        return item_iterator(new_item, parent.node);
    }
    
    
    /** 
     * @brief Insert a new tree node under the given node. 
     *
     * If the created node is the first child of `node`, then it will contain 
     * all of `node`'s `LeafItems`, otherwise it will be empty. (This guarantees 
     * that all leaf items have a parent node which is a leaf node).
     *
     * @return The inserted node.
     */
    node_iterator insert_child_node(const node_iterator& node, const NodeItem& obj) {
        NodeRef parent = node.node;
        NodeRef new_node;
        
        if (parent->n_children == 0) {
            NodeRef insert_pt = node_successor(parent);
            
            new_node = m_nodes.insert(insert_pt, Node(this, parent, obj));
            
            // take ownership of parent's items.
            new_node->items_first = parent->items_first;
            new_node->items_last  = parent->items_last;
            new_node->n_items     = parent->n_items;
            
            // update boundary.
            parent->child_first = parent->child_last = new_node;
        } else {
            // simple case: parent already contains node; 
            // insert point in straightforward (after the last child).
            NodeRef insert_pt = parent->child_last; ++insert_pt;
            new_node = m_nodes.insert(insert_pt, Node(this, parent, obj))
            
            // update last child.
            parent->child_last = new_node;
        }
        
        parent->n_children++;
        
        return new_node;
    }
    
    
    /**
     * @brief Insert a new tree node to the left of `insert_before`, under `parent`.
     *
     * If `insert_before` is the root node, or if `insert_before` is not a child (or end-child) of `parent`, then the tree is 
     * unaffected and `end()` is returned. Otherwise, the new node is returned.
     */
    node_iterator insert_sibling_node(
            const node_iterator& insert_before, 
            const node_iterator& parent, 
            const NodeItem& obj) {
        NodeRef n = insert_before.node;
        NodeRef p = parent.node;
        
        if (n == m_nodes.begin())                      return m_nodes.end();
        if (n->parent != p and parent.end().node != n) return m_nodes.end();
            
        NodeRef new_node = m_nodes.insert(n, Node(this, parent, obj));
        if (n == p->end()) {
            p->child_last = new_node;
        }
        if (n == p->begin()) {
            p->child_first = new_node;
        }
        p->n_children++;
        return new_node;
    }
    
    
    /**
     * @brief Insert new tree nodes to the left of `insert_before`, under `parent`.
     *
     * If `insert_before` is the root node, or if `insert_before` is not a child (or end-child) of `parent`, then the tree is 
     * unaffected and `end()` is returned. Otherwise, the new node is returned.
     *
     * @param insert_before A node_iterator pointing to the sibling just after the last new node to be inserted.
     * @param parent Parent node of the nodes to be created.
     * @param i_begin A forward iterator to the first NodeItem to be inserted.
     * @param i_end A forward iterator just beyond the last NodeItem to be inserted.
     */
    template <typename NodeItemIterator>
    node_iterator insert_sibling_nodes(
            const node_iterator& insert_before,
            const node_iterator& parent,
            const NodeItemIterator& i_begin,
            const NodeItemIterator& i_end) {
        NodeRef n = insert_before.node;
        NodeRef p = parent.node;
        
        if (n == m_nodes.begin())                      return m_nodes.end();
        if (n->parent != p and parent.end().node != n) return m_nodes.end();
            
        // todo: m_nodes.reserve(n_items);
        NodeRef new_node;
        index_t ct = 0;
        for (NodeItemIterator i = i_begin; i != i_end; ++i, ++ct) {
            new_node = m_nodes.insert(n, Node(this, parent, *i));
        }
        if (n == p->end()) {
            p->child_last = new_node;
        }
        if (n == p->begin()) {
            p->child_first = new_node;
        }
        p->n_children += ct;
        return new_node;
    }
    
    
    /**
     * @brief Split `node` into two sibling nodes.
     * 
     * A new node will be created to the left of `node`, and the items of `node` will be split 
     * between them according to the result of `compare(item, pivot)`: If the comparison is 
     * less than zero, the item will be moved to the left (new) node. Otherwise, it will
     * remain in the right (existing) node.
     *
     * @param node The node to split; a non-root node with no child nodes.
     * @param compare A callable object `compare(a,b)` which accepts a `LeafItem` as its left 
     *        argument and a `P` as its right argument, and returns a signed number 
     *        (negative for `a < b`, positive for `a > b`, zero for equality).
     * @param obj The internal node object to assign to the newly-created (low) node.
     *
     * @return The newly created sibling, or the off-end node if `node` is the root node or has child nodes.
     */
    template <typename Func, typename P>
    node_iterator split(node_iterator& node, Func compare, P pivot, NodeItem& obj) {
        if (node == m_nodes.begin() or node.nodes() > 0) return m_nodes.end();
        
        node_iterator new_node = insert_sibling_node(node, obj);
        NodeRef hi = node.node;
        NodeRef lo = new_node.node;
        
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
    }
    
    
    /**
     * @brief Remove an item from the tree which is a descendent of the node `ancestor`.
     *
     * If `item` is not in the subtree, the tree will be unchanged.
     *
     * Providing a close ancestor of `item` reduces the cost of finding `item`'s parent in the tree before removal.
     * In general, an `O(n * log(n))` search is performed to find the direct parent of `item`, where `n` is the
     * size of the subtree belonging to `ancestor`.
     *
     * @return An iterator pointing to the position of the item just beyond the one that was deleted. 
     * This will be `parent.items_end()` if the deleted item was the last one in its parent node.
     */
    item_iterator erase(const item_iterator& item, const node_iterator& ancestor) {
        NodeRef parent = find_parent(ancestor, item).node;
        
        if (parent != m_nodes.end()) return _erase(item.item, parent);
        
        return item_iterator(m_items.end(), m_nodes.begin());
        
        // todo: in key/value subclass, you can find `item` in log(n) time. 
        // be sure to override this fucker.
        
        // todo: no way to detect error condition of "item not found".
    }
    
    
    /**
     * @brief Remove an item from the tree.
     *
     * Performs an `n log(n)` search through the node from which `item` was obtained, where `n` is the number
     * of items in that node's subtree.
     *
     * @return An iterator pointing to the position of the item just beyond the one that was deleted, 
     * or `items_end()` if it was the last one in the tree.
     */
    item_iterator erase(const item_iterator& item) {
        return erase(item, item.parent);
    }
    
    
    /**
     * @brief Empty `node` of all child nodes and descendent items.
     */
    void clear(const node_iterator& node) {
        flatten(node);
        
        NodeRef n = node.node;
        
        if (n->n_items > 0) {
            // remove deleted items from ancestors.
            ItemRef del_first_item = p->items_first;
            ItemRef del_last_item  = p->items_last;
            ItemRef endpt = p->items_last; ++endpt;
            
            ItemRef next_item = m_items.erase(p->items_first, endpt);
            ItemRef prev_item = next_item; --prev_item;
            
            while (n != m_nodes.end()) {
                n->n_items -= n_deleted;
                // fix up the boundary items
                if (n->n_items == 0) {
                    n->items_first = n->items_last = m_items.end();
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
     * @brief Remove the subtree and its root at `node`, including all its leaf items.
     * 
     * `node` may not be the root node; if the root is given, the tree is unchanged.
     *
     * @return An iterator pointing to the position of the node following the deleted one (which may be `parent.end()`).
     */
    node_iterator erase(const node_iterator& node) {
        if (node.node->parent == m_nodes.end()) return false;
        
        clear(node);
        
        NodeRef x = node.node;
        
        NodeRef next_node = m_nodes.erase(x);
        NodeRef prev_node = next_node; --prev_node;
        
        // clean up the parents.
        NodeRef n = x->parent;
        
        // remove self from parent
        n->n_children--;
        if (n->n_children == 0) {
            // `x` was the only child
            n->child_first = n->child_last = m_nodes.end();
        } else if (n->child_first == x) {
            n->child_first = next_node;
        } else if (n->child_last == x) {
            n->child_last = prev_node;
        }
        
        return next_node;
    }
    
    
    /**
     * @brief Recursively remove all the children of `node`, and take ownership of all items beneath them.
     */
    void flatten(const node_iterator& node) {
        if (node->nodes() == 0) return; // nothing to do.
        NodeRef parent = node.node;
        
        // because of breadth-first order, this is exactly all the descendents of `parent`:
        m_nodes.erase(parent->child_first, node_successor(parent));
        
        parent->n_children = 0;
        parent->child_first = parent->child_last = m_nodes.end();
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
        if (not (a_n == m_nodes.end() and b_n == other.m_nodes.end()) return false;
        
        return true;
    }
    
    
    /// Tree inequality
    bool operator!=(const Tree<NodeItem, LeafItem>& other) const {
        return not (*this == other);
    }
    
    
protected:
    
    // erase when `parent` is known to be the exact parent of `item`; not just an ancestor.
    item_iterator _erase(const ItemRef& item, NodeRef parent) {
        ItemRef next_item = m_items.erase(item);
        ItemRef prev_item = next_item; --prev_item;
        
        // update the ancestors.
        while (parent != m_nodes.end()) {
            parent->n_items--;
            if (parent->n_items == 0) {
                parent->items_first = parent->items_last = m_items.end();
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
        
        return item_iterator(next_item, parent);
    }
    
    
    // use all the child counts to recompute the first/last 
    // references in internal nodes. useful when taking ownership
    // of another tree's internal node list.
    void recalculate_references() {
        NodeRef nod = m_nodes.begin();
        NodeRef begin_child = nod; ++begin_child;
        ItemRef begin_items = m_items.begin(); // next unassigned leaf item
        ItemRef last_item   = m_items.end();   // last / most recent assigned leaf item in the list
        
        for (; nod != m_nodes.end(); ++nod) {
            nod->child_first = nod->child_last = (nod->n_children > 0) ? begin_child : m_nodes.end();
            nod->items_first = nod->items_last = (nod->n_items    > 0) ? begin_items : m_items.end();
            
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
                if (n->n_items > 0) {
                    for (index_t i = 0; i < nod->n_items - 1; ++i) {
                        ++nod->items_last;
                    }
                    // move the item "free space"
                    begin_items = last_item = nod->items_last; ++begin_items;
                }
                    
                // fixup the items_end along the right edge of this subtree.
                NodeRef descendent = nod;
                for (NodeRef ancestor = nod->parent; ancestor != m_nodes.end(); descendent = ancestor, ancestor = ancestor->parent) {
                    if (ancestor->child_last == descendent) {
                        if (ancestor->n_items > 0) ancestor->items_last = last_item;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    
    
    // find the position of the node that follows `node` in list ("breadth-first") order.
    // the successor to `node` occurs just before the next subtree to the right.
    // an `O(k * log(n))` operation, on `n` the size of the tree, and `k` the arity.
    NodeRef node_successor(const NodeRef& node) const {
        NodeRef ancestor  = node->parent;
        NodeRef child     = node;
        NodeRef successor = m_nodes.end(); // if there is no subtree to the right, 
                                           // this is the conceptual successor.
        
        while (ancestor->parent != m_nodes.end()) {
            NodeRef end_child = ancestor->child_last; ++end_child;
            // look for a populated node to the right of us
            while (++child != end_child) {
                if (child->n_children > 0) {
                    successor = child->child_first;
                    break;
                }
            }
            if (successor != m_nodes.end()) {
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
    
    
    // find the first item that occurs after where the last item owned by the given node would go, in tree order.
    ItemRef item_successor(const NodeRef& node) const {
        for (NodeRef parent = node->parent; parent != m_nodes.end(); node = parent, parent = node->parent) {
            for (NodeRef end_c = parent->child_last; node != end_c; ++node) {
                if (node->n_items > 0) {
                    return node->items_first;
                }
            }
        }
        return m_items.end();
    }
    
}; // class Tree


/// @} // addtogroup storage


} // namespace geom
