#include <geomc/shape/shapedetail/IndexHelpers.h>
#include <geomc/Templates.h>

#include <cmath>
#include <list>

namespace geom {


// todo: std::list may be arbitrarily stupid with memory allocation. look into this.
// todo: owned std::list would allow a native, light swap() on items (aka splice?)
//       owned list would also permit reserve()
//       owned list could avoid copying unnecessary objects to newly constructed child nodes
//       > or consider going back to flat arrays. :[
//       > consider looking into how to make/use a std::allocator.
// todo: helper should first try to call dist2() on bounded objects, 
//       then revert to obj.bounds().dist2().
// todo: could make the bound a generic Concept type.
// todo: make a thin Container wrapper for i.getPoints() and i.rangeSearch()
// todo: range search
// todo: knn?
// todo: generic proximity search. some kind of distance<>(a,b), overlaps<a,b>, matches<a,b>, etc.
//       e.g. ray or segment intersect, generic shape intersect, disjoint tests, and so on.

// todo: come up with a portable alloca() and remove the arity restriction
//       you could do: fast_buffer<N>(n). always statically allocates N, 
//                     if n > N, malloc also and use that instead
// todo: xxx: fix clients of Rect to handle boundary issue
// todo: xxx: erase() on objects.
// todo: xxx: enforce that node_arity <= MAX_ARITY when inserting kids.
// todo: xxx: median pivot

/** @addtogroup shape
 *  @{
 */


/*****************************************
 * Enums                                 *
 *****************************************/

/// Strategy for choosing an axis along which to split a group of objects in a spatial index.
enum class KDAxisChoice {
    /// Divide a node along its longest axis.
    AXIS_LONGEST,
    /// Divide a node along the axis with the highest position variance.
    AXIS_HIGHEST_VARIANCE,
    /// Divide each node along an axis that cycles with tree depth, starting with the x (first) axis at the root node.
    AXIS_CYCLICAL
};


/// Strategy for choosing the pivot value when splitting a group of objects in a spatial index.
enum class KDPivotChoice {
    /// Split a node at the population average.
    PIVOT_MEAN,
    /// Split a node at the population median.
    PIVOT_MEDIAN
};


/// Strategy for insertion of an object into an existing tree.
enum class KDInsertionChoice {
    /// Insert objects into the node that will experience the smallest volume increase as a result of the addition.
    INSERT_SMALLEST_VOLUME_INCREASE,
    /// Insert objects into the node which is closest.
    INSERT_SHORTEST_DISTANCE
};


// if I had a cross-platform alloca(), I could instead obey the zero/one/infinity princicple. >:(
// I am going to use this space to gripe about c++ apologists who think "you can't do it, therefore you don't need it".
#define MAX_KDTREE_ARITY 16


/*****************************************
 * KDTree class                          *
 *****************************************/


/**
 * @brief A hierarchical spatial index.
 *
 * @tparam T Numeric type.
 * @tparam N Dimensionality of data.
 * @tparam Object Spatial object to be indexed. May be a `Vec<T,N>`, a `Bounded<T,N>`, or a 
 * `std::pair<K,V>` with `K` either of the former.
 * @tparam NodeData Optional data to store with each internal node of the tree.
 */ 
template <typename T, index_t N, 
          typename Object, 
          typename NodeData>
class KDTree {

private:
    
    /************************************
     * Internal type definitions        *
     ************************************/
    
    struct KDNode;
    
    typedef typename std::list<KDNode>::iterator KDNodeRef;
    typedef typename std::list<Object>::iterator KDDataRef;
    typedef detail::ShapeIndexHelper<T,N,Object> helper_t;
    
    struct KDNode {
        
        KDNodeRef  parent;
        KDNodeRef  child_begin;
        KDNodeRef  child_end;
        KDDataRef  objects_begin;
        KDDataRef  objects_end;
        index_t    nobjs;
        
        Rect<T,N>  bounds;
        NodeData   data;
    };

public:
    
    /**
     * An optionally-const iterator over the internal nodes of a KD tree.
     *
     * Dereferencing this iterator produces a `NodeData` object.
     */
    template <bool Const>
    class KDNodeIterator {
        
        friend class KDTree<T,N,Object,NodeData>;
        
        KDNodeRef node;
        
        KDNodeIterator(KDNodeRef n) : node(n) {}
        
        KDNodeIterator(const KDNodeIterator<false>& i) : node(i.node) {}
        
    public:
        
        typedef NodeData                                         value_type;
        typedef typename ConstType<NodeData,Const>::reference_t  reference;
        typedef typename ConstType<NodeData,Const>::pointer_t    pointer;
        typedef const NodeData&                                  const_reference;
        typedef KDNodeIterator<Const>                            iterator;
        typedef KDNodeIterator<Const>                            self_t;
        typedef typename ConstType<Rect<T,N>,Const>::reference_t bound_reference;
        typedef typename std::conditional<Const,
                            typename std::list<Object>::const_iterator,
                            typename std::list<Object>::iterator>::type   object_iterator;
        
        /// `+i`: Become first child
        inline self_t& operator+() {
            node = node->child_begin;
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
            KDNodeIterator tmp = *this;
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
        
        /// Get first child
        inline self_t begin() const {
            return node->child_begin;
        }
        
        /// Get last (off-end) child
        inline self_t end() const {
            return node->child_end;
        }
        
        /// Get first object inside this node
        inline object_iterator objects_begin() const {
            return node->objects_begin;
        }
        
        /// Get last (off-end) object in this node
        inline object_iterator objects_end() const {
            return node->objects_end;
        }
        
        /// Get bound
        inline bound_reference bound() {
            return node->bounds;
        }
        
    };
    
    
    /// Iterator over the internal tree nodes
    typedef KDNodeIterator<false> node_iterator;
    
    /// Const iterator over the internal tree nodes
    typedef KDNodeIterator<true>  const_node_iterator;
    
    /// Iterator over objects
    typedef typename std::list<Object>::iterator object_iterator;
    
    /// Const iterator over objects
    typedef typename std::list<Object>::const_iterator const_object_iterator;
    
    
    /// Structure encapsulating the tree balancing parameters
    struct KDStructureParams {
        /// Strategy for choosing an axis to split when subdividing a node
        KDAxisChoice axis;
        /// Strategy for choosing a pivot value when subdividing a node
        KDPivotChoice pivot;
        /// Strategy for insertion of an object into an existing node.
        KDInsertionChoice insert;
        /// Maximum number of children an internal node can have
        index_t node_arity;
        /// Maximum number of data objects a leaf node can have
        index_t leaf_arity;
    };
    
private:
    
    /************************************
     * Members                          *
     ************************************/
    
    // invariant after construction: root node exists.
    
    std::list<KDNode> nodes;
    std::list<Object> objects;
    KDStructureParams params;
    
    static const KDStructureParams DefaultParameters;
    
    // in general, whenever we update the structure of the tree,
    // we must be careful to:
    //   - update the relevant nodes' `nobjs`
    //   - update the relevant nodes' bounds
    //   - repair any invalidated iterators
    //   - rebalance the tree if a node has become too heavy or light
    
public:
    
    
    /************************************
     * Structors                        *
     ************************************/
    
    
    /// Construct an empty KDTree with the default structure parameters for this dimension.
    KDTree():params(DefaultParameters) {
        nodes.push_back(KDNode());
        KDNode&       n = nodes.front();
        n.parent        = nodes.end();
        n.child_begin   = nodes.end(); // no children
        n.child_end     = nodes.end();
        n.objects_begin = objects.begin();
        n.objects_end   = objects.end();
        n.nobjs         = 0;
    }
    
    
    /// Construct a KDTree initialized with `objs`. 
    KDTree(Object* objs, index_t nobjs, const KDStructureParams params=DefaultParameters) {
        //objects.reserve(nobjs);
        objects.insert(objects.begin(), objs, objs + nobjs);
        //nodes.reserve(std::ceil(std::log(nobjs, params.node_arity)));
        rebalance(params);
    }
    
    
    /// Construct a KDTree initialized with the objects contained in the interval [begin, end).
    template<typename ObjectIterator>
    KDTree(ObjectIterator begin, ObjectIterator end, const KDStructureParams params=DefaultParameters) {
        objects.insert(objects.begin(), begin, end);
        //nodes.reserve(std::ceil(std::log(objects.size(), params.node_arity)));
        rebalance(params);
    }
    
    
    /************************************
     * Functions                        *
     ************************************/
    
    
    /// Number of data objects in the tree.
    inline index_t nobjects() const {
        return nodes.begin()->nobjs;
    }
    
    
    /**
     * Return an iterator pointing at the root tree node. 
     * Use `+i` and `-i` to ascend and descend to the first-child and parent nodes
     * respectively. `++i` and `--i` navigate the tree in breadth-first order, and can
     * therefore be used to find the next and previous siblings.
     */
    inline node_iterator begin() {
        return nodes.begin();
    }
    
    
    /**
     * Return an iterator pointing beyond the last node in the tree; conceptually the parent of the root.
     */
    inline node_iterator end() {
        return nodes.end();
    }
    
    
    /**
     * Return a const iterator pointing at the root tree node. 
     * Use `+i` and `-i` to ascend and descend to the first-child and parent nodes
     * respectively. `++i` and `--i` navigate the tree in breadth-first order, and can
     * therefore be used to find the next and previous siblings.
     */
    inline const_node_iterator begin() const {
        return nodes.begin();
    }
    
    
    /**
     * Return a const iterator pointing beyond the last node in the tree in breadth-first order; 
     * conceptually the parent of the root.
     */
    inline const_node_iterator end() const {
        return nodes.end();
    }
    
    
    /**
     * Insert the given object into the tree under the given node.
     * The tree will be descended recursively until the best leaf is found.
     */
    object_iterator insert(const node_iterator& node, const Object& obj) {
        Rect<T,N> bnd = helper_t::bounds(obj);
        Vec<T,N>  ctr = helper_t::getPoint(obj);
        
        index_t depth;
        if (params.axis == KDAxisChoice::AXIS_CYCLICAL) {
            // compute depth
            depth = 0;
            for (KDNodeRef i = node.node->parent; i != nodes.end(); i = i->parent) {
                ++depth;
            }
        } else {
            depth = 1;
        }
        
        return insert_impl(node.node, obj, bnd, ctr, depth);
    }
    
    
    /**
     * Remove the object at the given iterator from the given leaf node.
     *
     * If the node is not a leaf or does not contain the object, no change
     * is made and the original iterator is returned. Otherwise, an iterator to the 
     * next object is returned.
     *
     * Both iterators may be invalidated by this operation.
     */
    object_iterator erase(const object_iterator& obj, const node_iterator& node) {
        KDNodeRef n = node.node;
        
        // verify node is a leaf
        if (n->child_begin != n->child_end) return obj;
        
        // verify `obj` is in this node
        bool found = false;
        Rect<T,N> bnd;
        for (KDDataRef i = n->objects_begin; i != n->objects_end; ++i) {
            if (i == obj) {
                found = true;
            } else {
                bnd |= helper_t::bounds(*i);
            }
        }
        
        if (not found) {
            return obj;
        } else {
            n->bounds = bnd;
            n->nobjs--;
        }
        
        // xxx todo: fixup anybody who used us as a first or a last.
        
        // save this, because we might invalidate n
        KDNodeRef parent = n->parent;
        
        // node is too empty. kill it.
        if (n->npts < std::max(2, params.leaf_arity / 4)) {
            // ? only one other sibling -> rebalance parent
            // ? node has multiple siblings -> reinsert??
            // xxx todo
        }
    }
    
    
    /**
     * Rebuild the given subtree according to its current balancing parameters, in `O(n log(n))` time.
     * If rebuilding from the root, it is more efficient to call `rebalance()` with no args.
     *
     * Any iterators pointing to descendents of the given node are invalidated.
     */
    void rebalance(const node_iterator& node) {
        flatten(node);
        index_t depth = -1;
        // for now, we only need the depth if we are using a cyclical axis scheme
        if (params.axis == KDAxisChoice::AXIS_CYCLICAL) {
            depth = 0;
            KDNodeRef n = node.node;
            for (KDNodeRef n = node.node; n != nodes.end(); n = n->parent) {
                depth++;
            }
        }
        buildTree(node.node, depth);
    }
    
    
    /**
     * Rebuild the entire tree.
     *
     * All iterators are invalidated.
     */
    void rebalance() {
        nodes.clear();
        nodes.push_back(KDNode());
        
        KDNode&       n = nodes.front();
        n.parent        = nodes.end();
        n.child_begin   = n.child_end = nodes.end(); // no children
        n.objects_begin = objects.begin();
        n.objects_end   = objects.end();
        n.nobjs         = objects.size();
        
        buildTree(nodes.begin(), 0);
    }
    
    
    /**
     * Change the balancing parameters of the tree, and rebuild it according to the new artings.
     */
    inline void rebalance(const KDStructureParams& newParams) {
        params = newParams;
        params.node_arity = std::min(params.node_arity, (index_t)MAX_KDTREE_ARITY);
        rebalance();
    }
    
    
    /**
     * Remove the given node, all of its children, and all the objects it contains.
     * Runs in `O(n)` time on the node count of the deleted subtree.
     *
     * The root node cannot be erased. If it is passed as an argument,
     * it is returned without deletion.
     *
     * Any iterators pointing to this node, any of its descendents, or the objects
     * contained therein are invalidated. If deleting this node leaves one remaining 
     * sibling, then the sibling will be absorbed into its parent and its iterator
     * will be invalidated.
     * 
     * @return The node after the deleted node, in breadth-first order.
     */
    node_iterator erase(const node_iterator& node) {
        KDNodeRef n = node.node;
        if (n == nodes.begin()) return node; // do not erase root
        
        clear(node); // empty this node
        
        KDNodeRef tail = n; ++tail;
        KDNodeRef parent = n->parent;
        
        // if we were the first child, parent needs a new favorite
        if (parent->child_begin == n) {
            parent->child_begin = tail;
        }
        
        // we were an uncle's last-child
        KDNodeRef fixup = n; --fixup;
        fixup = fixup->parent;
        if (fixup != nodes.end() and fixup->child_end == n) {
            fixup->child_end = tail;
        }
        
        tail = nodes.erase(node);
        
        KDNodeRef last_sibling = parent->child_end; --last_sibling;
        if (parent->child_begin == last_sibling) {
            // parent has one child. collapse it.
            bool invalid_tail = parent->child_begin == tail;
            KDNodeRef new_tail = flatten(parent);
            if (invalid_tail) tail = new_tail;
        }
        
        return tail;
    }
    
    
    /**
     * Empty the given node of all its objects, and delete all its child nodes.
     * Runs in `O(n)` time on the node count of the emptied subtree.
     *
     * Any iterators pointing to nodes or objects in the deleted subtree are invalidated.
     */
    void clear(const node_iterator& node) {
        KDNodeRef n = node.node;
        KDDataRef tail = n->objects_end;
        
        // ascend the tree, looking for nodes that refer to the deleted range
        // and close the hole we're about to make
        KDNodeRef fixup = n->parent;
        while (fixup != nodes.end() and fixup->objects_begin == n->objects_begin) {
            // ancestor's range began with deleted range
            if (fixup->objects_begin == fixup->objects_end) {
                fixup->objects_begin  = fixup->objects_end = tail;
            } else {
                fixup->objects_begin = tail;
            }
            // ancestor's sibling ended with deleted range
            KDNodeRef sibling = fixup; --sibling;
            if (sibling != nodes.end() and sibling->objects_end == n->objects_begin) {
                sibling->objects_end = tail;
                if (sibling->objects_begin = n->objects_begin) {
                    sibling->objects_begin = tail;
                }
            }
        }
        
        // the right edge of the previous sibling subtree ends with the deleted range.
        // fixup those nodes to point to the other side of the hole we just made.
        // (nothing to do if there is no subtree to the left of us)
        if (n->parent != nodes.end() and n != n->parent->child_begin) {
            fixup = n; --fixup; // previous sibling (we know from test above that it exists)
            while (true) {
                fixup->objects_end = tail;
                if (fixup->child_begin == fixup->child_end) break;
                fixup = fixup->child_end; --fixup;
            }
        }
        
        // remove child nodes
        flatten(node);
        
        // remove objects
        index_t removed_objs = n->nobjs;
        objects.erase(n->objects_begin, n->objects_end);
        n->objects_begin = n->objects_end = tail;
        n->bounds = Rect<T,N>();
        n->nobjs = 0;
        
        // update ancestors' bounds to reflect missing objects
        while (n->parent != nodes.end()) {
            n = n->parent;
            n->bounds = Rect<T,N>();
            for (KDNodeRef i = n->child_begin; i != n->child_end; ++i) {
                n->bounds |= i->bounds;
            }
            n->nobjs -= removed_objs;
        }
    }
    
    
    /** 
     * Delete all the child nodes of the given node, and assign all the 
     * objects below to it. Runs in `O(n)` time on the node count of the 
     * emptied subtree.
     *
     * Any iterators pointing to nodes in the deleted subtree are invalidated.
     */
    void flatten(const node_iterator& node) {
        KDNodeRef n = node.node;
        if (n->child_begin == n->child_end) return;
        for (KDNodeRef c = n->child_begin; c != n->child_end; c++) {
            flatten(c);
        }
        
        // we deleted a node that is some other node's child_end.
        // we need to suture up the hole that we made.
        KDNodeRef fixup = n.child_begin;
        fixup = (--fixup)->parent;
        KDNodeRef tail = nodes.erase(n.child_begin, n.child_end);
        if (fixup != nodes.end()) {
            if (fixup->child_begin == fixup->child_end) {
                fixup->child_begin  = fixup->child_end = tail;
            } else {
                fixup->child_end = tail;
            }
        }
        
        // we now have no children.
        n->child_begin = n->child_end = tail;
    }
    
    
    /** 
     * Insert a new tree node under the given node. If the created node is the
     * first child of `node`, then it will contain all of `node`'s objects,
     * otherwise it will be empty.
     *
     * @return The inserted node.
     */
    node_iterator insertChild(const node_iterator& node) {
        KDNodeRef n = node.node;
        KDNodeRef new_node;
        if (n->child_begin == n->child_end) {
            // creating the first child of `node`.
            new_node = nodes.insert(n->child_end, *n);
            new_node->child_begin = new_node->child_end = n->child_end;
            new_node->parent = n;
            n->child_begin = new_node;
            if (n->parent != nodes.end() and n->parent->child_end == n->child_end) {
                // are we creating a "first grandchild"?
                // this node is the first one not a direct child of node->parent.
                n->parent->child_end = new_node;
            }
        } else {
            // creating a new sibling
            new_node = nodes.insert(n->child_end, KDNode());
            new_node->parent        = n;
            new_node->child_begin   = new_node->child_end   = n->child_end;
            new_node->objects_begin = new_node->objects_end = n->objects_end;
        }
        return new_node;
    }
    
    
    /**
     * Return the object in the tree closest to the query point `q`.
     */
    const Object& nearest(const Vec<T,N>& p) const {
        KDDataRef ref;
        T big = std::numeric_limits<T>::max();
        nearestInNode(nodes.begin(), p, &ref, big, big);
        return *ref;
    }
    
    /**
     * Return the object in the tree closest to the query point `q`.
     */
    Object& nearest(const Vec<T,N>& p) {
        KDDataRef ref;
        T big = std::numeric_limits<T>::max();
        T best_d2    = big;
        T worst_best = big;
        nearestInNode(nodes.begin(), p, &ref, &best_d2, &worst_best);
        return *ref;
    }
    
    
    /// Get the parameters describing the current tree-balancing strategy.
    inline const KDStructureParams& getStructureParams() { return params; }


private:
    
    
    /**
     * Recalculate the bound of the given node by directly unioning the bounds
     * of all the objects inside.
     */
    Rect<T,N> recalculateBounds(KDNodeRef node) {
        Rect<T,N> bnd;
        for (KDDataRef i = node->objects_begin; i != node->objects_end; ++i) {
            bnd |= helper_t::bounds(*i);
        }
        node->bounds = bnd;
        return bnd;
    }
    
    
    /**
     * Make a new sibling of `split_child` (assumed a leaf) and divide its objects between them
     * according to the current balance parameters. Do not check if the node needs to be split,
     * or if splitting will exceed the desired node arity.
     */
    void splitNode(KDNodeRef split_child, index_t depth) {
        Vec<T,N> mean;
        Vec<T,N> dim;
        Vec<T,N> var;
        
        index_t split_axis = 0;
        T pivot;
        
        // choose a split axis
        switch (params.axis) {
            case KDAxisChoice::AXIS_LONGEST:
                dim = split_child->bounds.getDimensions();
                goto CHOOSE_BIGGEST_AXIS;
            case KDAxisChoice::AXIS_CYCLICAL:
                split_axis = depth % N;
                break;
            case KDAxisChoice::AXIS_HIGHEST_VARIANCE: {
                for (KDDataRef i = split_child->objects_begin; i != split_child->objects_end; ++i) {
                    mean += helper_t::getPoint(*i);
                }
                for (KDDataRef i = split_child->objects_begin; i != split_child->objects_end; ++i) {
                    Vec<T,N> p = helper_t::getPoint(*i) - mean;
                    var += p * p;
                }
                var /= std::max<index_t>(1, split_child->nobjs - 1); // variance with Bessel's correction
                dim = var;
                
            CHOOSE_BIGGEST_AXIS:
                split_axis = 0;
                T big = dim[0];
                for (index_t i = 1; i < N; i++) {
                    if (dim[i] > big) {
                        big = dim[i];
                        split_axis = i;
                    }
                }    
                break;
            }
        }
        
        // choose a pivot
        switch (params.pivot) {
            case KDPivotChoice::PIVOT_MEAN:
                if (params.axis == KDAxisChoice::AXIS_HIGHEST_VARIANCE) {
                    // already calculated mean above
                    pivot = mean[split_axis];
                } else {
                    // calculate the mean now
                    for (KDDataRef i = split_child->objects_begin; i != split_child->objects_end; ++i) {
                        pivot += helper_t::getPoint(*i)[split_axis];
                    }
                    pivot /= split_child->nobjs;
                }
                break;
            case KDPivotChoice::PIVOT_MEDIAN:
                // xxx: todo: this
                break;
        }
        
        // make a new child with no objects and no children
        KDNodeRef new_sibling = insertChild(split_child->parent).node;
        
        // partition the objects; create bbox
        split_child->bounds = Rect<T,N>();
        index_t nobjs = split_child->nobjs;
        KDDataRef middle = new_sibling->objects_end = split_child->objects_end;
        for (KDDataRef i = split_child->objects_begin; i != middle; ) {
            Vec<T,N> v = helper_t::getPoint(*i);
            if (v[split_axis] > pivot) {
                // move this object to the new (upper) sibling
                std::swap(*i, *--middle);
                split_child->nobjs--;
                new_sibling->bounds |= helper_t::bounds(*i);
                // `i` now contains an unchecked object; don't increment.
            } else {
                split_child->bounds |= helper_t::bounds(*i);
                ++i;
            }
        }
        split_child->objects_end = new_sibling->objects_begin = middle;
        new_sibling->nobjs = nobjs - split_child->nobjs;
    }
    
    
    /**
     * Recursively sort the objects in `node` into subtrees according to the
     * current balancing parameters in `O(n log(n))` time.
     */
    void buildTree(KDNodeRef node, index_t depth) {
        // root node needs to calc its box
        if (depth == 0) {
            recalculateBounds(node);
        }
        
        if (node->nobjs > params.leaf_arity) {
            
            KDNodeRef split_child = insertChild(node).node;
            index_t n_added = 1;
            
            while (n_added < params.node_arity and split_child != node->child_end) {
                // continue splitting the oldest node 
                // until we've either reached node arity or run out of splittable nodes
                if (split_child->nobjs > params.leaf_arity) {
                    splitNode(split_child, depth);
                    ++n_added;
                }
                ++split_child;
            }
            
            // recurse
            for (KDNodeRef i = node->child_begin; i != node->child_end; ++i) {
                buildTree(i, depth + 1);
            }
            
        }
    }
    
    
    /**
     * Insert the given object into `node`, descending recursively until a leaf is found.
     */
    object_iterator insert_impl(
            KDNodeRef node, 
            const Object& obj, 
            const Rect<T,N>& bnd,
            const Vec<T,N>& ctr,
            index_t depth) {
        
        node->bounds |= bnd;
        node->nobjs++;
        
        if (node->child_begin != node->child_end) {
            KDNodeRef best;
            T best_metric = std::numeric_limits<T>::max();
            for (KDNodeRef i = node->child_begin; i != node->child_end; ++i) {
                T metric;
                switch (params.insert) {
                    case KDInsertionChoice::INSERT_SHORTEST_DISTANCE:
                        metric = i->bounds.dist2(ctr);
                        break;
                    case KDInsertionChoice::INSERT_SMALLEST_VOLUME_INCREASE:
                        T orig_vol = i->bounds.volume();
                        T new_vol  = (i->bounds | bnd).volume();
                        metric = new_vol - orig_vol;
                        break;
                }
                if (metric == 0) {
                    metric = -1 / i->bounds.getCenter().dist2(ctr);
                }
                if (metric < best_metric) {
                    best = i;
                    best_metric = metric;
                }
            }
            object_iterator out = insert_impl(best, obj, bnd, ctr);
        } else {
            objects.insert(node->objects_end, obj);
            
            // subdivide if necessary
            if (node->nobjs > params.leaf_arity) {
                index_t parent_arity = 0;
                KDNodeRef parent = node->parent;
                if (parent != nodes.end()) {
                    // count parent's children
                    for (KDNodeRef i = parent->child_begin; i != parent->child_end; ++i) {
                        ++parent_arity;
                    }
                    // new sibling, or new subtree?
                    if (parent_arity >= params.node_arity) {
                        buildTree(node, depth);
                    } else {
                        splitNode(node, depth);
                    }
                } else {
                    // root node cannot have new silbings
                    buildTree(node, depth);
                }
            }
        }
    }
    
    
    /**
     * Find the object nearest `p` in the given node.
     */
    void nearestInNode(KDNodeRef node, const Vec<T,N>& p, KDDataRef* nearest, T* best_d2, T* worst_best) const {
        struct SearchInfo {
            KDNodeRef n;
            T best_case;
        };
        index_t child_ct   = 0;
        index_t recurse_ct = 0;
        
        // todo: With a portable alloca(), this could be arbitrarily large or small.
        // this is the only line of code where leaf arity is actually restricted. 
        // but instead we are forced to use the maximum sane amount of memory in all cases.
        // (ironically, the argument against alloca is that it "uses too much stack").
        SearchInfo buf[MAX_KDTREE_ARITY];
        
        for (KDNodeRef i = node->child_begin; i != node->child_end and recurse_ct <= MAX_KDTREE_ARITY; ++i, ++child_ct) {
            T box_d2 = i->bounds.dist2(p);
            T box_worst_d2 = detail::worst_case_nearest2(p, i->bounds);
            if (box_d2 < *best_d2 and box_d2 <= *worst_best) {
                buf[recurse_ct++] = {i, box_d2};
            }
            *worst_best = std::min(*worst_best, box_worst_d2);
        }
        // recurse
        while (--recurse_ct >= 0) {
            // subsequently added nodes may have better worst cases than we knew about before,
            // or subsequent recursions may have improved the best estimate, rendering our
            // checks obsolete and allowing us to discard more nodes. check again plz.
            if (buf[recurse_ct].best_case < *best_d2 and buf[recurse_ct].best_case <= *worst_best) {
                nearestInNode(buf[recurse_ct].n, p, nearest, best_d2, worst_best);
            }
        }
        // is leaf?
        // seach objects for actual nearest obj.
        if (child_ct == 0) {
            for (KDDataRef i = node->objects_begin; i != node->objects_end; ++i) {
                T d2 = helper_t::dist2(*i, p);
                if (d2 < *best_d2) {
                    *best_d2 = d2;
                    *nearest = i;
                }
            }
        }
    }
    
    
}; // KDTree class


/// @} // addtogroup shape


/*****************************************
 * Static members                        *
 *****************************************/


template <typename T, index_t N, typename Object, typename NodeData>
const typename KDTree<T,N,Object,NodeData>::KDStructureParams KDTree<T,N,Object,NodeData>::DefaultParameters = 
{
    KDAxisChoice::AXIS_LONGEST,
    KDPivotChoice::PIVOT_MEAN,
    KDInsertionChoice::INSERT_SMALLEST_VOLUME_INCREASE,
    std::min(2 << N, MAX_KDTREE_ARITY),        // internal node arity
    std::min(2 << N, MAX_KDTREE_ARITY)         // leaf node arity
};


} // namespace geom