#include <geomc/geomc_defs.h>
#include <geomc/Templates.h>
#include <list>

// todo: permit key/value structure.
//       nodes represent an abstract kind of set membership. node.contains(key). 
//       (maybe call it boundary so people don't retardedly nest data structures).
// todo: after key/value structure, permit a distance metric; allow KNN.
// todo: consider factoring the above simply as free functions which implement algorithms.
// todo: specialize for when leaf items or node items are void
// todo: smarter/more optimized underlying data structure
// todo: "flat array" variant tree. abstract NodeRef and ItemRef to be ptrs.
// todo: figure out a way to put node_successor, aka subtree_last into node iterator?
// todo: implement find(item_test(LeafItem) -> bool, node_test(NodeItem)=lambda x:true)
//       a generic algorithm for walking the tree to search for an item.
//       node_test() returns `true` if the node might possibly contain the item.
//       (templating over the functions would allow functors and inlining)
// todo: might be good to support re-rooting a tree at a child
// todo: might be good to support subtree-adoption from another tree


// xxx: coredump in node_successor
//      > made a change; need to test.
//      > change was ascent check against ::end()

// xxx: note that deleting a node invalidates all item iterators
//      beneath it, in that they might have considered the deleted
//      node to be their "subtree root". this seems undesireable,
//      especially because it would not be the case if the subtree root
//      were not stored.

namespace geom {
    
/** @addtogroup storage
 *  @{
 */

/**
 * @brief A dynamic tree of arbitrary arity.
 *
 * Each node in the tree may store an associated `NodeItem` object. 
 * A node may have any number of child nodes.
 * 
 * Leaf nodes may own a list of `LeafItem` objects. `LeafItem`s are kept
 * such that all the `LeafItems` in a subtree are contiguous, and internal nodes
 * provide iterators to the first and last `LeafItem`s in their respective subtrees.
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
        
        index_t  n_items;    // i.e. all descendents
        index_t  n_children; // i.e. direct children
        
        NodeItem data;
    };
    
    
    // internally, nodes are kept in "sibling-first" order:
    //   - visit a node
    //   - visit its next sibling, if there is one
    //   - visit its first child
    // as such, all direct children of a node are contiguous, and come 
    // after (though not necessarily *immediately* after) their parent.
    // all the descendents of a leftward child come before all the 
    // descendents of the child to its right. Therefore, any given subtree
    // (excluding its root) is contiguous. It is this recursive pattern:
    
    // [node] ... [c1 c2 c3 c4 ...] 
    //    [children of c1] [descendents of c1...] 
    //    [children of c2] [descendents of c2...] 
    //    ...
      
    std::list<Node>     nodes;
    
    // items are stored in depth-first order. 
    // all items in a subtree are contiguous.
    
    std::list<LeafItem> items;
    
public:
    
    template <bool Const>
    class ItemIterator {
        
        friend class Tree<NodeItem,LeafItem>;
        
        typedef ItemIterator<Const> self_t;
        // reference to value type
        typedef typename ConstType<LeafItem,Const>::reference_t reference;
        // pointer to value type
        typedef typename ConstType<LeafItem,Const>::pointer_t   pointer;
        
        ItemRef item;
        NodeRef parent; // an ancestor of the item. may or may not be "close" 
                        // depending on how the item was found. this can be used 
                        // to accelerate tree search/mutation actions on an 
                        // item iterator.
        
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
    
    public: // xxx WHY
        
        friend class Tree<NodeItem,LeafItem>;
        
        NodeRef node;
        
        NodeIterator(NodeRef n) : node(n) {}
        
         // can always promote/cast to const:
        NodeIterator(const NodeIterator<false>& i) : node(i.node) {}
        
    public:
        
        // container concept typedefs:
        
        /// The tree's `NodeItem` type.
        typedef NodeItem                                         value_type;
        /// A (possibly const) reference to a `NodeItem`.
        typedef typename ConstType<NodeItem,Const>::reference_t  reference;
        /// A (possibly const) pointer to a `NodeItem`.
        typedef typename ConstType<NodeItem,Const>::pointer_t    pointer;
        /// A const reference to a `NodeItem`.
        typedef const NodeItem&                                  const_reference;
        /// Iterator over child nodes; same as self type.
        typedef NodeIterator<Const>                              iterator;
        // self type
        typedef NodeIterator<Const>                              self_t;
        
        /// Iterator over `LeafItem`s.
        typedef ItemIterator<Const> item_iterator; // iterator for leaf items
        
        
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
         * @brief Get last (off-end) object in this subtree. 
         * It is invalid to increment or dereference this iterator.
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
        Node&       root = m_nodes.front();
        root.parent      = m_nodes.end();
        root.child_first = m_nodes.end(); // no children
        root.child_last  = m_nodes.end();
        root.items_first = m_items.end(); // no items
        root.items_last  = m_items.end();
        root.n_items     = 0;
        root.n_children  = 0;
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
     * @brief Return a const iterator pointing beyond the last node in the tree.
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
     * @brief Return the first node in the subtree rooted at `root` 
     * (excluding `root` itself).
     */
    inline node_iterator subtree_begin(const node_iterator& root) {
        return root->child_first;
    }
    
    
    /**
     * @brief Return the last (off-end) node in the subtree rooted at `root`.
     */
    inline node_iterator subtree_end(const node_iterator& root) {
        return node_successor(root.node);
    }
    
    
    /**
     * @brief Return the first node in the subtree rooted at `root` 
     * (excluding `root` itself).
     */
    inline const_node_iterator subtree_begin(const const_node_iterator& root) const {
        return root->child_first;
    }
    
    
    /**
     * @brief Return the last (off-end) node in the subtree rooted at `root`.
     */
    inline const_node_iterator subtree_end(const const_node_iterator& root) const {
        return node_successor(root.node);
    }
    
    
    /**
     * @brief Exhaustively search for the direct parent node of item `i` in 
     * the subtree rooted at `ancestor`.
     *
     * If `i` is not in the subtree under `node`, then return `end()`.
     */`
    node_iterator find_parent(
            const const_node_iterator& ancestor, 
            const const_item_iterator& i) const {
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
     * @brief Insert an item into the given node.
     *
     * The new item will be inserted before the item at `insert_before`.
     * `insert_before` must belong to the given node, otherwise the behavior 
     * is undefined. It is permissible for `insert_before` to be the node's 
     * `items_begin()` or `items_end()`.
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
     * @param obj New leaf item to be inserted.
     *
     * @return An iterator to the first newly placed item if any were placed;
     * `items_end()` otherwise.
     */
    item_iterator insert(
            const node_iterator& parent,
            const item_iterator& insert_before,
            const LeafItem& obj) {
        // todo: it would be great if we could protect against the user
        //   providing an insert_pt that doesn't belong to the parent
        
        // must insert to a leaf node
        if (parent.node_count() > 0) {
            return items_end();
        }
        
        // retrieve parent and insert point
        NodeRef node = parent.node;
        ItemRef insert_pt = insert_before.item;
        if (node.n_items == 0) {
            insert_pt = item_successor(node);
        }
        
        // insert the item
        ItemRef new_item = m_items.insert(insert_pt, obj);
        
        // ascend the tree, updating the item ranges of the ancestors
        ItemRef precursor_item = new_item; --precursor_item;
        while (node != m_nodes.end()) {
            if (node->items_first == insert_pt or node->n_items == 0) {
                node->items_first = new_item;
            }
            if (node->items_last == precursor_item or node->n_items == 0) {
                node->items_last = new_item;
            }
            node->n_items += 1;
            node = node->parent;
        }
        
        return new_item;
    }
    
    
    /**
     * @brief Insert multiple leaf items into the given node.
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
    item_iterator insert(
            const node_iterator& parent,
            const item_iterator& insert_before,
            const LeafItemIterator  begin_item,
            const LeafItemIterator& end_item,
            index_t* new_item_count=nullptr) {
        // todo: it would be great if we could protect against the user
        //   providing an insert_pt that doesn't belong to the parent
        
        // must insert to a leaf node
        if (parent.node_count() > 0) {
            return items_end();
        }
        if (begin_item == end_item) {
            // nothing to be done
            return items_end();
        }
        
        // retrieve parent and insert point
        NodeRef node = parent.node;
        ItemRef insert_pt = insert_before.item;
        if (node.n_items == 0) {
            insert_pt = item_successor(node);
        }
        
        // insert the first item
        ItemRef first_new_item = m_items.insert(insert_pt, *begin_item);
        ItemRef last_new_item  = first_new_item;
        ++begin_item;
        
        // insert the remaining items
        index_t item_ct = 1;
        for (LeafItemIterator i = begin_item; i != end_item; ++i) {
            last_new_item = m_items.insert(insert_pt, *i);
            ++item_ct;
        }
        
        // ascend the tree, updating the item ranges of the ancestors
        ItemRef precursor_item = first_new_item; --precursor_item;
        while (node != m_nodes.end()) {
            if (node->items_first == insert_pt or node->n_items == 0) {
                node->items_first = first_new_item;
            }
            if (node->items_last == precursor_item or node->n_items == 0) {
                node->items_last = last_new_item;
            }
            node->n_items += item_ct;
            node = node->parent;
        }
        
        if (new_item_count) *new_item_count = item_ct;
        return first_new_item;
    }
    
    /**
     * @brief Insert a new tree node to the left of `insert_before`, 
     * under `parent`.
     *
     * If the parent node was previously empty, then its new first child will
     * contain all of the parent's `LeafItem`s (this is to ensure those leaf
     * items all have a parent which is a leaf node). Otherwise, the new node
     * will be empty of any `LeafItem`s.
     *
     * If `insert_before` is the root node, or if `insert_before` is not 
     * a child (or end-child) of `parent`, then the tree is unaffected 
     * and `end()` is returned. Otherwise, the new node is returned.
     *
     * Invalidates all `end()` `node_iterator`s.
     *
     * @param parent Parent node of the new node to be created.
     * @param insert_before A `node_iterator` pointing to the sibling just
     * after the new node.
     * @param obj The `NodeItem` to be associated with the new node.
     * @return The newly created node, or `end()` if the node could not be 
     * created.
     */
    node_iterator insert_child_node(
            const node_iterator& parent,
            const node_iterator& insert_before,
            const NodeItem& obj) {
        NodeRef p = parent.node;
        NodeRef n = insert_before.node;
        bool parent_empty = p->n_children == 0;
        
        // check invalid cases
        if (n == m_nodes.begin())                      return m_nodes.end();
        if (n->parent != p and parent.end().node != n) return m_nodes.end();
        
        // set the insert point for empty nodes
        if (parent_empty) {
            n = node_successor(p)
        }
        
        // create the new node
        NodeRef new_node = m_nodes.insert(n, Node(this, parent, obj));
        
        // take ownership of the parent's items, if we are its first child
        // (every item must belong to a leaf node, and the parent
        // has just become an internal node).
        if (parent_empty) {
            new_node->items_first = p->items_first;
            new_node->items_last  = p->items_last;
            new_node->n_items     = p->n_items;
        }
        
        // adjust child endpoints if necessary
        if (insert_before.node == parent.end().node) {
            p->child_last = new_node;
        }
        if (insert_before.node == p.begin().node) {
            p->child_first = new_node;
        }
        p->n_children++;
        return new_node;
    }
    
    
    /**
     * @brief Insert new tree nodes to the left of `insert_before`, 
     * under `parent`. Return the first newly created node.
     *
     * If the parent node was previously empty, then its new first child will
     * contain all of the parent's `LeafItem`s (this is to ensure those leaf
     * items all have a parent which is a leaf node). All other new nodes
     * will be empty of any `LeafItem`s.
     *
     * If `insert_before` is the root node, or if `insert_before` is not 
     * a child (or end-child) of `parent`, then the tree is unaffected and 
     * `end()` is returned. Likewise, if no new nodes were inserted, 
     * `end()` is returned. Otherwise, the first new node is returned. 
     * 
     * Invalidates all `end()` `node_iterator`s.
     *
     * @param parent Parent node of the new nodes to be created.
     * @param insert_before A node_iterator pointing to the sibling just 
     * after the last new node to be inserted.
     * @param i_begin A forward iterator to the first NodeItem to be inserted.
     * @param i_end A forward iterator just beyond the last NodeItem to be inserted.
     * @return The first newly created node, or `end()` if no new nodes were created.
     */
    template <typename NodeItemIterator>
    node_iterator insert_child_nodes(
            const node_iterator& parent,
            const node_iterator& insert_before,
                  NodeItemIterator  i_begin,
            const NodeItemIterator& i_end) {
        NodeRef p = parent.node;
        NodeRef n = insert_before.node;
        bool parent_empty = p->n_children == 0;
        
        // check invalid cases
        if (n == m_nodes.begin())                      return m_nodes.end();
        if (n->parent != p and parent.end().node != n) return m_nodes.end();
        if (i_begin == i_end)                          return m_nodes.end();
        
        // set the insert point for empty nodes
        if (parent_empty) {
            n = node_successor(p);
        }
        
        // create the first node, and keep a reference to it
        // todo: m_nodes.reserve(n_items);
        NodeRef first_new_node = m_nodes.insert(n, Node(this, parent, *(i_begin++)));
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
            last_new_node = m_nodes.insert(n, Node(this, parent, *i));
        }
        
        // adjust child endpoints if necessary
        if (insert_before.node == parent.end().node) {
            p->child_last = last_new_node;
        }
        if (insert_before.node == parent.begin().node) {
            p->child_first = first_new_node;
        }
        p->n_children += ct;
        return first_new_node;
    }
    
    
    /**
     * @brief Split `node` into two sibling nodes.
     * 
     * A new node will be created to the left of `node` under the same parent, 
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
     * @param node The node to split; a non-root node with no child nodes.
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
    node_iterator split(node_iterator& node, Func compare, P pivot, NodeItem& obj) {
        if (node == m_nodes.begin() or node.node_count() > 0) return m_nodes.end();
        
        node_iterator new_node = insert_child_node(node.node->parent, node, obj);
        NodeRef hi = node.node;
        NodeRef lo = new_node.node;
        
        // xxx: can the formulation of this be flipped so that
        //      the new node is higher (rightward / more recent),
        //      rather than lower? it should act more like an `append`.
        
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
     * @brief Remove an item from the subtree rooted at `ancestor`.
     *
     * If `item` is not in the subtree, the tree will be unchanged.
     *
     * Providing a close ancestor of `item` reduces the cost of finding 
     * `item`'s parent in the tree before removal. In general, an `O(n * log(n))` 
     * search is performed to find the direct parent of `item`, where `n` is the
     * size of the subtree belonging to `ancestor`.
     *
     * Invalidates the iterator to this item, as well as any `end()`
     * `item_iterator`s.
     *
     * @return An iterator pointing to the position of the item just beyond 
     * the one that was deleted. This will be `parent.items_end()` if the 
     * deleted item was the last one in its parent node.
     */
    item_iterator erase(
            const item_iterator& item, 
            const node_iterator& ancestor) {
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
     * Performs an `n log(n)` search through the subtree from which `item` 
     * was obtained, where `n` is the number of items in that subtree.
     *
     * Invalidates the iterator to this item, as well as any `end()`
     * `item_iterator`s.
     *
     * @return An iterator pointing to the position of the item just 
     * beyond the one that was deleted, or `items_end()` if it was the 
     * last one in the tree.
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
            ItemRef del_first_item = n->items_first;
            ItemRef del_last_item  = n->items_last;
            ItemRef endpt = n->items_last; ++endpt;
            
            // ItemRef next_item = m_items.erase(n->items_first, endpt);
            // ItemRef prev_item = next_item; --prev_item;
            
            // xxx: check list::erase spec to make sure below is equivalent
            //      (i.e. return value position)
            ItemRef next_item;
            index_t n_deleted;
            for (next_item = n->items_first; next_item != endpt; ++next_item) {
                m_items.erase(next_item);
                ++n_deleted;
            }
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
     * @return An iterator pointing to the position of the node following the 
     * deleted one (which may be `parent.end()`).
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
     * @brief Recursively remove all the children of `node`, and take 
     * ownership of all items beneath them.
     */
    void flatten(const node_iterator& node) {
        if (node->node_count() == 0) return; // nothing to do.
        NodeRef parent = node.node;
        
        // because of breadth-first order, 
        // this is exactly all the descendents of `parent`:
        m_nodes.erase(parent->child_first, 
            (parent));
        
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
        if (not (a_n == m_nodes.end() and b_n == other.m_nodes.end())) return false;
        
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
        
        // in essence we "build" the tree by simply taking ownership of 
        // the nodes/items that already exist. we walk over the nodes in 
        // tree order, which means that child nodes are taken in order too.
        // this spares us recursion (i.e. no possibility of stack-blowing) 
        // and some logic keeping track of descent/ascent.
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
                        ancestor != m_nodes.end(); 
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
    
    
    // find the position of the node that follows `node` in list ("sibling-first") order.
    // the successor to `node` occurs just before the next subtree to the right.
    // an `O(k * log(n))` operation, on `n` the size of the tree, and `k` the arity.
    // this is different from `++node` because `++node` may descend into the 
    // current node's subtree; we want to exactly hop over it.
    NodeRef node_successor(const NodeRef& node) const {
        NodeRef ancestor  = node->parent;
        NodeRef child     = node;
        NodeRef successor = nodes_end().node; // if there is no subtree to the right,
                                              // this is the conceptual successor.
        
        while (ancestor != m_nodes.end()) {
            NodeRef off_end_child = ancestor->child_last; ++off_end_child;
            // look for a populated node to the right of us
            while (++child != off_end_child) {
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
    
    
    // find the first item that occurs after where the last item owned by 
    // the given node would go, in tree order.
    ItemRef item_successor(NodeRef node) const {
        for (NodeRef parent = node->parent; 
                parent != m_nodes.end(); 
                node = parent, parent = node->parent) {
            NodeRef end_c = parent->child_last; ++end_c;
            for (; node != end_c; ++node) {
                if (node->n_items > 0) {
                    return node->items_first;
                }
            }
        }
        return items_end().item;
    }
    
}; // class Tree


/// @} // addtogroup storage


} // namespace geom
