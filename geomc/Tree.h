#include <geomc/Templates.h>

// todo: function(s) to move an item/child from one node to another?
//       (may require custom list)
// todo: permit key/value structure.
//       nodes represent an abstract kind of set membership. node.contains(key). 
//       (maybe call it boundary so people don't retardedly nest data structures).
// todo: after key/value structure, permit a distance metric; allow KNN.
// todo: specialize for when leaf items or node items are void?


// todo: return next iterator on erase(). can be compared to node.end().

namespace geom {
    
/** @addtogroup storage
 *  @{
 */

template <typename NodeItem, typename LeafItem>
class Tree {
    
    struct Node;
    
    typedef typename std::list<Node>::iterator     NodeRef;
    typedef typename std::list<LeafItem>::iterator DataRef;
    
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
        DataRef  items_first;
        DataRef  items_last;
        
        index_t  n_items;
        index_t  n_children;
        
        NodeItem data;
    };
    
    /// Iterator over items
    typedef typename std::list<LeafItem>::iterator item_iterator;
    
    /// Const iterator over items
    typedef typename std::list<LeafItem>::const_iterator const_item_iterator;
    
    std::list<Node>     nodes;
    std::list<LeafItem> items;
    
public:
    
    /**
     * An optionally-const iterator over the internal nodes of a tree.
     *
     * Dereferencing this iterator produces a `NodeItem` object.
     */
    template <bool Const>
    class NodeIterator {
        
        friend class Tree<NodeItem,LeafItem>;
        
        NodeRef node;
        
        NodeIterator(NodeRef n) : node(n) {}
        
        NodeIterator(const NodeIterator<false>& i) : node(i.node) {}
        
    public:
        
        typedef NodeItem                                         value_type;
        typedef typename ConstType<NodeItem,Const>::reference_t  reference;       // reference to value type
        typedef typename ConstType<NodeItem,Const>::pointer_t    pointer;         // pointer to value type
        typedef const NodeItem&                                  const_reference; // const ref to value type
        typedef NodeIterator<Const>                              iterator;        // self type (container concept)
        typedef NodeIterator<Const>                              self_t;          // self type
        typedef typename std::conditional<Const,
                    typename std::list<Object>::const_iterator,
                    typename std::list<Object>::iterator>::type  item_iterator; // iterator for leaf items
        
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
        inline index_t nodes() const {
            return node->n_children;
        }
        
        /// Number of leaf items in this subtree
        inline index_t items() const {
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
        
        /// Get first object inside this node
        inline item_iterator items_begin() const {
            return node->items_first;
        }
        
        /// Get last (off-end) object in this node
        inline item_iterator items_end() const {
            item_iterator tmp = node->items_last;
            return ++tmp;
        }
        
    };
    
    
    /// Iterator over the internal tree nodes
    typedef NodeIterator<false> node_iterator;
    
    /// Const iterator over the internal tree nodes
    typedef NodeIterator<true>  const_node_iterator;
    
    
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
     * Construct a tree having a single root node, and copy the
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
    
    
    /// Number of leaf items in the tree.
    inline index_t items() const {
        return m_nodes.begin()->n_items;
    }
    
    
    /// Number of internal nodes in the tree.
    inline index_t nodes() const {
        return m_nodes.begin()->n_children + 1; // include root, of course
    }
    
    
    /**
     * Return an iterator pointing at the root tree node. 
     * Use `+i` and `-i` to ascend and descend to the first-child and parent nodes
     * respectively. `++i` and `--i` navigate the tree in breadth-first order, and can
     * therefore be used to find the next and previous siblings.
     */
    inline node_iterator begin() {
        return m_nodes.begin();
    }
    
    
    /**
     * Return an iterator pointing beyond the last node in the tree; conceptually the parent of the root.
     */
    inline node_iterator end() {
        return m_nodes.end();
    }
    
    
    /**
     * Return a const iterator pointing at the root tree node. 
     * Use `+i` and `-i` to ascend and descend to the first-child and parent nodes
     * respectively. `++i` and `--i` navigate the tree in breadth-first order, and can
     * therefore be used to find the next and previous siblings.
     */
    inline const_node_iterator begin() const {
        return m_nodes.begin();
    }
    
    
    /**
     * Return a const iterator pointing beyond the last node in the tree in breadth-first order; 
     * conceptually the parent of the root.
     */
    inline const_node_iterator end() const {
        return m_nodes.end();
    }
    
    
    /**
     * Return parent node of item `i`, if `i` belongs to the subtree at `n`,
     * by performing a depth-first search in `n * log(n)` time.
     *
     * If `i` is not in the subtree `n`, then return `end()`.
     */
    NodeRef find_parent(const NodeRef& n, const item_iterator& i) const {
        // todo: do not recurse. make a loop, to avoid possible stack overflow.
        if (n->n_children == 0) {
            item_iterator end_item = n->items_last; ++end_item;
            for (item_iterator j = n->items_first; j != end_item; ++j) {
                if (j == i) return n;
            }
        } else {
            NodeRef end_child = n->child_last; ++end_child;
            for (NodeRef child = n->child_first; child != end_child; ++child) {
                NodeRef found = find_parent(child, j);
                if (found != m_nodes.end()) return found;
            }
        }
        return m_nodes.end();
    }
    
    
    /**
     * Insert the given item into the tree under the given node.
     * 
     * `node` must be a leaf node; i.e. one with no child nodes, so that placement of `obj` is not ambiguous.
     * If `node` is not a leaf node, then the tree is unchanged and the off-end item_iterator is returned instead.
     */
    item_iterator insert(const node_iterator& parent, const LeafItem& obj) { // todo: move semantics?
        if (parent.nodes() > 0) return items.end();
        
        // this is a bit complex and special case-y, but fundamentally it is O(log(n)).
        // we will walk the tree vertically from the insertion node at least once and at most twice.
        // most of this complexity has to do with handling the case where the new item becomes
        // the new first_item or last_item of its ancestor(s). furthermore, the bulk of this is only
        // paid when the first item is added to an empty node.
        
        item_iterator new_item;
        NodeRef ancestor = parent.node;
        
        if (parent.items() > 0) {
            // basic/common case: node already has items.
            item_iterator old_item = ancestor->items_first;
            item_iterator insert_pos = old_item; ++insert_pos;
            new_item = m_items.insert(insert_pos, obj);
            
            // walk up the tree updating the boundary (if necessary) and the item count.
            while (ancestor != m_nodes.end()) {
                ancestor->n_items++;
                if (ancestor->items_last == old_item) {
                    ancestor->items_last = new_item;
                }
                ancestor = ancestor->parent;
            }
        } else {
            // adding the first item to an empty node.
            // find the first ancestor with childen
            NodeRef cur_node = parent.node;
            while (cur_node->n_items == 0 and cur_node->parent != m_nodes.end()) {
                ancestor = cur_node;
                cur_node = cur_node->parent;
            }
            
            // find a populated subtree adjacent to the inserted item's ancestor.
            // to ensure that items are contiguous in a tree, we must place the new item between them.
            NodeRef populated_left  = m_nodes.end();
            NodeRef populated_right = m_nodes.end();
            bool hit_child = false;
            // walk the children of cur_node looking for a populated subtree.
            for (NodeRef sibling = cur_node->child_first; true; ++sibling) {
                if (sibling == ancestor) {
                    hit_child = true; // child abuse
                    continue;
                }
                if (sibling->n_items > 0) {
                    // found a populated subtree.
                    if (not hit_child) {
                        // it's on the left side of the ancestor.
                        populated_left = sibling;
                    } else {
                        // it's on the right side of the ancestor.
                        populated_right = sibling;
                        break; // can't discover any new info now.
                    }
                }
                if (sibling == cur_node->child_last) break;
            }
            
            // place the new item between a populated subtree and the current one.
            bool empty_root = false;
            if (populated_left != m_nodes.end()) {
                // populated subtree is to the left of us; place new item to the right of it
                item_iterator tmp = populated_left->child_last; ++tmp;
                new_item = m_items.insert(tmp, obj);
            } else {
                if (populated_right != m_nodes.end()) {
                    // populated subtree is to the right of us; place new item to the left of it
                    new_item = m_items.insert(populated_right->child_first, obj);
                } else {
                    // empty root node.
                    new_item = m_items.insert(m_items.end(), obj);
                    cur_node->items_first = cur_node->items_last = new_item;
                    empty_root = true;
                }
            }
            
            // update the boundaries of the parent subtree, if necessary.
            // we may have inserted an item on the "outside" of the parent's first/last child.
            // (because we have at least one sibling to the left or right of us with items,
            // this can only happen on one of the two boundaries).
            if ((populated_left == m_nodes.end() or populated_right == m_nodes.end()) and not empty_root) {
                bool left_edge = populated_left == m_nodes.end(); // was it the left edge that was empty?
                item_iterator old_edge = (left_edge) ? (cur_node->items_first) : (cur_node->items_last);
                while (cur_node != m_nodes.end()) {
                    NodeRef& update = (left_edge) ? (cur_node->items_first) : (cur_node->items_last);
                    if (update == old_edge) {
                        update = new_item;
                    } else {
                        break;
                    }
                    cur_node = cur_node->parent;
                }
            }
            
            // walk upwards from the item's parent.
            // - update the item counts
            // - set the boundary items for all newly-nonempty ancestors of the placed item
            ancestor = parent.node;
            while (ancestor != m_nodes.end()) {
                ancestor->n_items++;
                if (ancestor->items_first == m_items.end()) {
                    ancestor->items_first = ancestor->items_last = new_item;
                }
                ancestor = ancestor->parent;
            }
        }
        
        return new_item;
    }
    
    
    /** 
     * Insert a new tree node under the given node. If the created node is the
     * first child of `node`, then it will contain all of `node`'s objects,
     * otherwise it will be empty. (This guarantees that all leaf items have a 
     * parent node which is a leaf node).
     *
     * @return The inserted node.
     */
    node_iterator insert_node(const node_iterator& node, NodeItem&& obj) {
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
            // insert point in straighforward (after the last child).
            NodeRef insert_pt = parent->child_last; ++insert_pt;
            new_node = m_nodes.insert(insert_pt, Node(this, parent, obj))
            
            // update last child.
            parent->child_last = new_node;
        }
        
        parent->n_children++;
        
        return new_node;
    }
    
    
    /**
     * Remove an item from the tree which belongs to the node `ancestor`.
     * If `item` is not in the subtree, the tree will be unchanged.
     *
     * Providing a close ancestor of `item` reduces the cost of finding `item`'s parent in the tree before removal.
     */
    bool erase(const item_iterator& item, const node_iterator& ancestor) {
        NodeRef parent = find_parent(ancestor.node);
        
        if (parent != m_nodes.end()) _erase(item, parent);
        
        // todo: in key/value subclass, you can find `item` in log(n) time. 
        // be sure to override this fucker.
    }
    
    
    /**
     * Remove an item from the tree.
     *
     * Performs an O(n) search through the tree to find `item`'s parent.
     */
    bool erase(const item_iterator& item) {
        return erase(item, this->begin());
    }
    
    
    /**
     * Remove all child nodes and all descendent items from `node`.
     */
    void clear(const node_iterator& node) {
        flatten(node);
        
        NodeRef n = node.node;
        
        if (n->n_items > 0) {
            // remove deleted items from ancestors.
            item_iterator del_first_item = p->items_first;
            item_iterator del_last_item  = p->items_last;
            item_iterator endpt = p->items_last; ++endpt;
            
            item_iterator next_item = m_items.erase(p->items_first, endpt);
            item_iterator prev_item = next_item; --prev_item;
            
            while (n != m_nodes.end()) {
                n->n_items -= n_deleted;
                if (n->items_first == del_first_item) {
                    n->items_first = next_item;
                }
                if (n->items_last == del_last_item) {
                    n->items_last = prev_item;
                }
                n = n->parent;
            }
        }
    }
    
    
    /**
     * Remove the subtree with root `node`, including all the leaf items.
     * 
     * `node` may not be the root node; if the root is given, the tree is unchanged.
     */
    bool erase(const node_iterator& node) {
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
    }
    
    
    /**
     * Recursively remove all the children of `node`, and take ownership of all items beneath them.
     */
    void flatten(const node_iterator& node) {
        if (node->nodes() == 0) return; // nothing to do.
        NodeRef parent = node.node;
        
        // because of breadth-first order, this is exactly all the descendents of `parent`:
        m_nodes.erase(parent->child_first, node_successor(parent));
        
        parent->n_children = 0;
        parent->child_first = parent->child_last = m_nodes.end();
    }
    
    
protected:
    
    // erase when `parent` is known to be the exact parent of `item`; not just an ancestor.
    void _erase(const item_iterator& item, NodeRef parent) {
        item_iterator next_item = m_items.erase(item);
        item_iterator prev_item = next_item; --prev_item;
        
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
    }
    
    
    // find the position of the node that follows `node` in list (breadth-first) order.
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
    
};


/// @} // addtogroup storage


} // namespace geom
