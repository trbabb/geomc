#include <geomc/Tree.h>
#include <geomc/CircularBuffer.h>

// todo: should we make BoundingHierarchy a base class, and then say that
//       user implementations are partial specializations that inherit from that?
//       - this would allow us to template the partial ordering type of insertion_quality
//       - it would also allow the bound(), combine(), etc. functions to be statically inlined
//       - however, the class would have to be curiously recurring and that is not elegant for the user.
// todo: how shall the user make direct edits to the items? is that allowed?
//       if the edit affects the bounding, it will break things.
//       but if we want to make an edit that *doesn't* affect the bounding, how do we do it?
//       --> store pointers.
//       --> ...although this still has the same gotcha.

// XXX: We want to separate the notion of whether a tree is editable and
//      whether its contents (items, internal nodes) are editable.
//      We want to let the user fuck with internal node data, including
//      by walking the tree. but if the iterators are const... not allowed.
//      how to fix...?

// todo: create BoundingScheme<Item,Boundary>, make specializations for various combinations.
//       have BoundingHierarchy take a template param for a bounding scheme which defaults
//       to BoundingScheme. users can provide a different class which doesn't specialize
//       or subclass boundingscheme if they want.
//       also, have bounding scheme provide combine(item, bound) which defaults to 
//       combine(bound(item), bound). also if type(item) == type(bound) then a specialization
//       needs to be made which doesn't (ambiguously) overload combine() this way.

namespace geom {

template <typename Item, typename Boundary, typename RegionData=void>
class BoundingHierarchy {
    
    struct Node {
        Boundary    bounds;
        RegionData  data;
    };
    
public:
        
    typedef Tree<Node, Item>                                tree_t;
    typedef ConstSubtree<Node,Item>                         const_node_iterator;
    typedef Subtree<Node,Item>                              node_iterator;
    typedef typename ConstSubtree<Node,Item>::item_iterator const_item_iterator;
    typedef typename Subtree<Node,Item>::item_iterator      node_iterator;
    
protected:
    
    tree_t  _tree;
    index_t _node_arity;
    index_t _item_arity;
    
public:
    
    BoundingHierarchy():
        _node_arity(8),
        _item_arity(8) {}
    
    /// Return a well-fitting Boundary that completely contains Item `i`.
    virtual Boundary bound(const Item& i) = 0;
    
    /// Return a Boundary which fully contains both `b0` and `b1`.
    virtual Boundary combine(const Boundary& b0, const Boundary& b1) = 0;
    
    /// Return `true` if Item `i` lies on or within Boundary `b`.
    virtual bool contains(const Boundary& b, const Item& i) = 0;
    
    /**
     * @brief Return a score representing how well Boundary `b` might accommodate the addition of Item `i`.
     * The insertion algorithm will try to minimize this score when inserting new items.
     *
     * For example, a scoring scheme might return the increase in volume to `b` after adding `i`.
     */
    virtual double insertion_quality(const Boundary& b, const Item& i) = 0;
    
    // todo: ^ would be better to pass a node_iterator so we can count how many things are inside?
    //       also: can we template the return type? all we need is a thing with partial ordering.
    
    /**
     * @brief Called when a node is created.
     *
     * When the tree is rebuilt, this function is called on child nodes
     * before it is called on parents.
     */
    virtual void node_created(const_node_iterator node) {}
    
    /// Called when a node's child count changes.
    virtual void node_changed(const_node_iterator node) {}
    
    /// Called when an item is inserted into the hierarchy.
    virtual void item_inserted(const_node_iterator external_node, const_item_iterator inserted_item) {}
    
    /// Called when an item is removed from the hierarchy.
    virtual void item_removed(const_node_iterator external_node, Item removed_item, const_item_iterator next_position) {}
    
    
    const_item_iterator insert(const Item& i) {
        // xxx todo
    }
    
    const_item_iterator erase(const_item_iterator i) {
        // xxx todo
    }
    
    const_node_iterator erase(const_node_iterator n) {
        // xxx todo
    }
    
    const_item_iterator find(const Item& i) const {
        // xxx todo
    }
    
    void rebuild(const_node_iterator i) {
        // remove constness. the user is not allowed to mutate our tree, but we are:
        node_iterator n = _tree.subtree(i);
        
        // we're redoing everything from the ground up:
        if (n.node_count() > 0) {
            n.flatten();
        }
        
        // xxx: todo: handle case where split produces an empty node, and the
        //      sibling produces an infinite loop.
        // xxx: todo: put the sibling of the split node back on the pile
        
        // split ourselves
        if (n.item_count() > _item_arity) {
            // make a new child identical to ourselves, put it on the stack
            CircularBuffer<node_iterator, 16> buf;
            buf.push_back(n.insert_child_node(*n);)
            
            // split into siblings until we are full, or all the siblings are small enough.
            while (buf.size() > 0 and n.node_count() < _node_arity) {
                node_iterator child = buf.pop_front();
                if (child.item_count() > _item_arity) {
                    // xxx how to choose a pivot?
                    //     beware that bounds and node data are not computed at this time.
                    //     note that using median/quickselect and only a comparator
                    //     obliterates the need for pivot-picking...
                    //     will also guarantee against "degenerate split" case
                    //        however: how would the caller choose an axis across which to split...?
                    node_iterator new_node = child.split(); // xxx fill in
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
                for (auto child = n.begin(); child != n.end(); ++c) {
                    rebuild(child);
                    rebound_this_level(child);
                    node_created(child);
                }
            } else {
                // split failed to create any new nodes.
                // remove the redundant node.
                // do not try to split again, it will also fail.
                n.erase(n.begin());
            }
        }
    }
    
    
    const tree_t& tree() const {
        return _tree;
    }
    
    // query_iterator query(Item i) const
    
protected:
    
    void rebound_this_level(node_iterator n) {
        if (n.node_count() > 0) {
            // use the bounds of the immediate children, assumed already correct.
            n->bounds = bound(*i.begin());
            for (auto i = ++begin(); i != n.end(); ++i) {
                n->bounds = combine(n->bounds, bound(*i));
            }
        } else if (n.item_count() > 0) {
            // leaf node. compute the bounds from the items.
            // (note: we should never have an empty node, but we checked for safety).
            n->bounds = bound(*n.items_begin());
            for (auto i = ++n.items_begin(); i != n.items_end(); ++i) {
                n->bounds = combine(n->bounds, bound(*i));
            }
        }
    }
    
    
}; // end BoundingHierarchy

} // end namespace geom