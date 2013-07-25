/*
 * DynamicExpr.h
 *
 *  Created on: Aug 16, 2011
 *      Author: tbabb
 */

//fun things that this could generalize:
//  image processing
//  audio processing
//  REYES
//  array processing/visualization, i.e. matlab

//competitors to investigate:
//  labview, simulink

/*
 * The current problem: execution context; evaluation order (topo sort)
 * 
 * Questions:
 *   - PRIMARY QUESTIONS:
 *     x how do we handle while(true)?
 *       - recursion
 *       - recursion-wrapping loop construct
 *         - a function described/drawn like a film frame (with sprocket holes?) that can connect 
 *           to the previous/next iteration
 *           - a switch that chooses between the output of "next N iterations" and direct output
 *             - this is the normal recursive way to do it
 *           - a switch that controls whether data is pushed to the next iteration, or allows
 *             direct output to flow past (foreach cell).
 *             - this kind of violates the directionality of data flow elsewhere: it's "pushing"
 *               vs. "reaching". 
 *               - that might be OK. You want the metaphor to convey what's happening.
 *             - what if we want some outputs now, but need to keep iterating for others?
 *             - how does this map to recursion?
 *               - almost directly. make the loop into a function. the final outputs are the fn outputs.
 *                 simply construct an if node under the "next iteration" representation, for each output.
 *         - a "next iteration" object which is placed at the bottom of every function?
 *       - function-building by key:value mapping
 *     x how do we handle axes of mismatching shape? e.g. (1, 2, 5) <--???--> (7, dynamic, 4, 22)
 *       > that would raise an error, because axis 0 lengths 1 and 7 can't be correlated.
 *         > message: Cannot join inputs from '%s' and '%s' in axis %d: '%s' has %d items; '%s' has %d items
 *                    Cannot join inputs from '+' and 'func5' in axis 0: '+' has 1 item; 'func5' has 7 items
 *       > if it were (5,2,6) <--> (5,dynamic,6,7) then axis 0 would correlate, axis 1 would be checked at run time,
 *         axis 2 would correlate, and axis 3 would see axes 0-2 as constant.
 *         > a user could also insert a new constant axis as any depth
 *   - TYPES:
 *     - three types of input: typed, function, dynamic. they should be drawn differently. (e.g. function
 *       connections are fatter). 
 *       - what about dimensionality?
 *       - what about constness? (per dimension)?
 *       - what about dimension dynamic-ness?
 *   - LIGHTNESS:
 *     - these should be as light as possible, suitable for constructing many, rapidly, at runtime.
 *   - COMPATIBILITY WITH EXPR: will connect()ing DynExprs result in the same structure as connect()ing static Exprs?
 *     - how can this be for constant expressions, which have input type void?
 *     - should ConstExprs lie about their input type, and say it is void when in fact it is something else?
 *     - should ConstExprs have specialized eval() functions, which we explicitly check for at the evaluation stage?
 *       - that check sounds slow.
 *     - important point: it should still be possible to use exprs as simple, lightweight functions
 *   - GRID DERIVATIVES: what about derivatives? "Filter regions" might be important for image processing, e.g.
 *     - requires an understanding of the dimensionality of the input.
 *     - should be built-in to language, or should a region be a type to operate on?
 *     - filtering type? (sinc? cubic? blackman? rect?)
 *     - should be possible to build this pretty easily. dimensionality is built in.
 *   - INPUT NAMES: are inputs indexed by name or position? both?
 *     - we'll say position for now, as naming causes some scoping issues
 *   - UNFILLED INPUTS: how are unfilled inputs handled?
 *     - Expressions are created using mutable prototypes which can't be directly used
 *     - The prototype is then flattened into a usable, immutable expr
 *     - Exprs/prototypes carry a complete, complex syntax tree, analagous to a function definition.
 *   - INPUT FINDING: how are inputs represented?
 *     - perhaps a linked list to unfilled leaf nodes?
 *   - EVALUATION ORDER: lazy evaluation-- sometimes it makes sense to evaluate up the tree, rather than down.
 *     e.g. conditionals; convolutions. How is this handled?
 *     - "up vs. down" is not the right way to think about it. Really we're talking about which
 *       nodes need to be evaluated before which. Conditionals, compositions, and selections mean that internal nodes
 *       take precedence over leaves.
 *       - it's not leaves vs. internal nodes; it's left vs. right branch.
 *     - compositions and convolutions are special in that they determine the domain over which to operate.
 *       So we need the concept of "dynamic domains". 
 *       - Is this really special? How is it different from the normal generation of domains/flow of 
 *         data?
 *       - How to implement it: 
 *         - input domain is evaluated first (i.e., space of the final image)
 *         - new domain is evaluated as a function of input domain
 *         - sampled function is evaluated over computed domain.
 *         - computed domain might be shorter or longer than the input domain.
 *           - longer for a convolution which goes off the edge, e.g.
 *           x shorter for a conditional, where some inputs are false.
 *             - "shorter" domains are those which are constant over 1+ axes
 *           x each bucket [cell] will have a domain of length >= 1. These domains
 *             should be collapsed and flattened.
 *             - BUT: mapping of generated buckets to their source buckets must be maintained (?)
 *             - this is a form of memoizing
 *           - domain-collapsing should happen in the case of a convolution, e.g.
 *             to avoid double-evaluation of the same sample point.
 *             - "collapsing" == repeat sample points are removed and only computed once
 *             - in some algorithms (separable convolutions?) this collapse might happen
 *               implicitly, in closed-form. will we ever waste time collapsing domains
 *               that an explicit (c++) implementation would avoid with more careful coding?
 *             - we do not have truly intelligent knowledge of which samples are going to be re-used.
 *               we might still re-evaluate per-grid. special-pupose code might be smart enough to 
 *               know when to "cache" data (in a var, e.g.) or discard it.
 *             > perhaps it might make sense for the user to do this explicitly; or give them a node to use
 *           - I like the idea of placing an evaluation cache manually.
 *           - the problem is how to arbitrarily re-combine longer domains back into one cell. 
 *             (reduce? iterated apply() of a function? these aren't terribly intutitive)
 *             - what about dynamic-programming type problems, where each output depends on C_n-1, C_n-2, ... ?
 *             - consider that the "distance" may be totally random. i.e. we might need item C_n-(some_dynamic_expr)
 *             - we might even need values across buckets
 *               - this is getting ugly again
 *               - what happened to "cached" values?
 *               - what about creating a node which turns a continuous function into a grid cache which can be sampled?
 *               - what about the IDisplace case?
 *           - really all you need to do in order to create new sub-domains is specify a length, per cell 
 *             (possibly constant) of the new sub-array. from there it can be transformed using standard
 *             techniques to the proper type of domain (from indecies). Will have to understand the concept
 *             of constancy across sub-domains as distinct from constancy across the cells which own those
 *             sub-domains, as well as the dimensionality of the sub-domain.
 *          - SUB-ARRAYS AS DIMENSIONS
 *            - it makes sense to think of having a sub-array the same as creating a new function
 *              with n+1 axis inputs. 
 *            - we then define ways to collapse axes. i.e by summation, reduction, iterated function application, etc.
 *            - i.e. for convolution, (x,y) -> f becomes (x, y, x', y') -> f
 *            - sub arrays are basically the only way to create inner loops, so this better f'in make sense and be
 *              intuitive and powerful.
 *              - except recursion, but that's not intuitive for many people 
 *   - CULLING:
 *     x notion of cullif(dataexpr, boolexpr)
 *       - deprecated in favor of a smarter if()
 *       - drops data buckets where true
 *     - keeping domains/data correlated is tricky
 *       - this is where "all domains are int" becomes interesting/handy
 *       - keep a list of indexes that are "live", and compact the data.
 *       - buckets that are not live in both inputs are culled
 *       - keep a "background" expr that is constant; not re-evaluated/iterated over by each new node. then used to fill in the gaps at the end
 *       - or perhaps we just need an if() that splits the buckets if one side is uniform and the other isn't
 *         - I kind of like this the best, it seems cleaner.
 *         - should probably do some checks, like number of const items. if is < 2, no reason to split.
 *         - evaluation order not obvious. 
 *           - Must evaluate the conditional expression before I know which buckets to cull
 *           - may not yield a culling at all if both sides are varying
 *           - which conditionals in the tree take priority?
 *   - CORRELATED DOMAINS: it seems we can't enforce correlated domains. This would break the ability to abstract expressions
 *     with inputs that will be filled later.
 *     - should we autocorrleate them when all leaves match domains, or should we make the user explicitly do it?
 *     * really, the question is: If we are requiring correlated domains, how do we expose multiple parameters?
 *     - we are only enforcing correlated domains in that in order to evaluate an expr, there must be a single input which
 *       maps to the integers.
 *   x USAGE AND SAMPLING: this raster-based evaluation method. let's think about it in the context of image processing.
 *     - you are thinking of having pixel values live in an unstructured "bag". How do we get X,Y?
 *       We may need that for convolutions, or displacement, etc. Fundamentally an image is a function
 *       of (x,y,c) -> val. we'd be destroying that?
 *       - No. color transformations are functions of vec3 -> vec3 (color to color), not vec2 -> vec3 (position to color).
 *         Thus the composition of an image with a color transformation is: vec2 -> vec3 -> vec3, which remains
 *         (vec2 -> vec3), an image.
 *       - How would a convolution work, though?
 *       - could this still be done if we model an image as vec3->val or vec4->val?
 *         - yes, via nodes that slice input?
 *   - EVALUATION OWNERSHIP: evaluation engine which walks the tree; manages caches?
 *     - "context" pointer instead, from which you can request caches; settings (i.e. whether to autocorrelate)?
 *       - might be okay, except the context would have to guess about the evaluation order, since it doesn't
 *         own the process. that would be important for figuring out whether a cache will be dirty/available/needed
 *         in the future.
 *   x NON-VARYING PARAMETERS: Nodes will need uniform-only parameters in (hopefully rare) cases; see PerlinNoise's seed param, which makes
 *     no sense to vary spatially.
 *     - but might make sense to vary over entries in a buffer, e.g. "per frame"
 *     - should be up to the user to decide
 *     - something has to tell the system to create new instances. i.e. the Expr subclass would have to declare 
 *       all its "one-per-instance" parameters, or make the copies internally (less desirable)
 *       - perlin is a node of input int, and an output of type Expr<Vec<T,N> -> T>
 *       - this can be wrapped so that we see a two-arg function of Vec<T,N>, int
 *   - PERFORMANCE: undesirable speed hits:
 *     - virtual indirection between DynamicExpr and Expr
 *     - type-check on well-formed expr trees
 *     - function call overhead of DynamicExpr->eval()
 *     - function call overhead of Expr->_impl_eval() for fnptr wrappers
 *     x multiple evaluation of uniform expressions when we've broken the execution into
 *       multiple buckets
 *       - going to make a cache of these
 *       - beware of stuff like recursive calls
 *   x COMPOSITION: with old-style exprs, composition is modeled as a node. should this be continued, or should
 *     nodes be inherently aware of their inputs?
 *     - NO: We need to distinguish between a definition and a usage. Usages occur in many places and need to
 *           know where their input data comes from. definitions occur once and are used in many places.
 *     - note that composition is currently a binop, with inputs F and G.
 *     - note that there is more than one way to think of this:
 *       - a function compose(), which takes two functions and returns a third function
 *       - a function compose(), which takes two functions and an input, and returns data. 
 *         - this one seems more stupider
 *     - the graphical tree would thus reflect the expression structure,
 *       not necessarily the data flow.
 *       - note: you could choose to draw composition nodes differently.
 *       - how much of an intuitive difference does this make?
 *       - which is better?
 *         - Old way (as a node):
 *           - Nodes are not aware of their inputs
 *           - Connecting an expr to another input does not alter it (safer; no side effects)
 *           - undesirable: non-unique representations for E = (A -> B -> C -> D)
 *             - e.g., (A -> B) -> (C -> D) or A -> ((B -> C) -> D); etc.
 *             - Thus when you ask for the parent of E, you might get (A -> B) if E is (C -> D); or 
 *               (A -> B -> C) if E is D. 
 *               - this could probably be accounted for in the parent-finding function
 *                 - New, unique exprs would probably have to be constructed, however.
 *               - an underlying and meaningful structure would be hidden by the graphical model.
 *               - Alternatively, a mostly unimportant distinction is introduced into graphical interaction.
 *         - New way (nodes keep pointers to data sources): 
 *           - Less modular: Different objects could not share a sub-procedure; data flow is fixed and non-variable all the way to
 *             the source.
 *             - Exprs would be altered as they are passed around 
 *             - Duplicates would have to be made in order for exprs to be shared and connected in different ways.
 *               - Thus connecting an input to your function costs O(number of nodes in the function), and
 *                 that's just stupid.
 *             - This is mostly fine for graphical interaction, but annoying for code interaction.
 *             - The most problematic issue is that changing/adding an input to a node
 *               could cause it to change its input type underneath you!
 *             - The in-code difficulties make me wary about the flexibility of the graphical implementation anyway.
 *               o Be specific?
 *           - Tree re-structuring has side-effects for all holders of a sub-expression
 *             - This is the case anyway with expression simplification, unless the process creates a duplicate.
 *           - Shake does it this way; each node is given a pointer to its inputs.
 *         - Consider: Exprs are "things to do". They are static. Perhaps the tree should be a list? Like RPN?
 *           - Trees and lists are equivalent. 
 *           - There is also an element of permutation: Certain branches are independent of each other and can
 *             begin execution simultaneously. (RPN encodes a specific order).
 *           - The most important thing to encode is the flow of information. From there, we can
 *             determine which procedures to invoke, and in what order. 
 *           - It should be easy to snip out branches and sub-graphs of this flow chart.
 *           - The graph should be lightweight.
 *           - The graph should be modular.
 *           - The graph is an acyclic digraph.
 *   - AXES: 
 *     - need to understand (and possibly distinguish between) input space of the function, number of raster axes,
 *       and dynamically expanded sub-domains.
 *       - is there a way we can map all of these concepts onto each other?
 *     - your notion of a static vector may not be the best way of dealing with multiple dimensions, since their
 *       dimensionality is fixed at compile time.
 *     - I like the idea of being able to change which axis is "major" 
 *       - does that make sense?
 *         - yes, it determines which axes we evaluate first. like:
 *           - calc the same pixel for every frame, then do the next pixel?
 *           - calc the whole image for one frame, then do the next frame?
 *             - no, that's not quite right. doesn't matter which axis is "major" inside a block;
 *               the entire block gets evaluated at once.
 *             - if you wanted to see a whole chunk of frame before the next one, you'd set the block
 *               size in the "time" axis to 1. A value of 5 would evaluate frames in chunks of 5. e.g.
 *               maybe 16x16x5 chunks of pixels.
 *             - so really what "major axis" means is the order of *block* evaluation.
 *               * so you could say "do all x first" or "do all frame 1 first". but what if you 
 *                 wanted to do a smaert 'out of order' eval, like your clever renderman idea?
 *         - should be totally possible to reorder these.
 *         - "dimensionality" might be the wrong presentation. We'd be building in the concept of
 *            holding one input the same while another varies.
 *            - it's important to be able to see what we're varying without having to
 *              choose/assign axis numbers.
 *            - maybe that should be done for you.
 *            - how does this map to loops?
 *     * how is constancy saved across blocks?
 *       - allocate buffer of len sufficient to hold necessary unique values
 *       - release value when the last block to use it is done
 *       - mem required will depend on # blocks per axis and size of other axes
 *     * what does block size mean in the context of varying length (ragged) axes?
 *       - you will be able to tell that an axis is ragged when the len param that feeds
 *         into its vary() node is not const across the parent axis.
 *     - named axes would be nice
 *       - ability to rename/remap
 *       - problem: global namespaces-- an outer function might choose the same axis name as an inner function
 *     - a table of axes | names | constancy
 *       - may not know "outer" axes, if, for example, our function is being called from an outer loop.
 *         - therefore might not know there's a dimension mismatch until we call the function
 *         - user might have to declare dimensionality
 *         - might require uniform dimensionality of all inputs? does this reduce functionality?
 *         * that's bad. it also means we have to deeply check dimension matching every time we wire up a function
 *     - lets say that each vector component becomes an input. then:
 *       - axis-constancy is nicely handled
 *       - cross-component operations get ugly.
 *       - example: an image becomes a function of (x,y,channel,time). 
 *         - you can partially fill each of those args to get a function of fewer dimensions
 *         - e.g. f(1,1,0,?) is a function describing how the red channel of pixel(1,1) changes over time.
 *         -      f(?,?,?,5) is an image (frame 5) which can be sampled.
 *         -      f(1,1,?,5) is a pixel at position 1,1 on frame 5
 *       - should be easy and intuitive to break one or more axes into a sub-problem
 *       - should be possible to sample buffers
 *       - this all comes back to operating on longer ranges
 *         - the ONLY reason dimensionality matters is:
 *           - what's in an inner loop vs. outer loop
 *           - what's constant for a loop vs. what's not.
 *         - there is a difference between just having two inputs vs. having two axes
 *           - sure, the two inputs can vary independently, but:
 *           - having two axes implies that one input is held constant while another varies, and
 *             vice versa.
 *         - notion of expandAxis, which holds its inputs constant per cell, generates a list of i, and accepts a function
 *           of the other inputs + i to fill the new cells
 *           - perhaps call it "array"? or "loop"? "vary"?
 *           - perhaps modeled as a node. it takes the dimensionality of its input and increases it by 1.
 *             - every joining input which does not pass through it is treated as constant across that axis.
 *             - it returns an integer which is the index (in the new dimension) of the current cell.
 *           - must be possible to index across cells
 *             - beware the dumb IDisplace case. Do you want to look at cells, or sample a function?
 *             - perhaps a node which remaps its input data as an ND function to be sampled? with output dimensionality
 *               of the sample point data?
 *               - but what of f(x,y) -> (a,b) vs. f(x,y,i) -> b? i.e. think about how you would represent a function of
 *                 (x,y) -> (a,b) and sample it in this manner.
 *               - how do we deal with sample inputs of differing dimensionality? 
 *           - dimensionality of data must be known (i.e. either constant or dynamic, i.e an expr) at every node
 *           - maybe it would be helpful to think of things this way:
 *             - all we are doing is keep track of regions of 1d data that are uniform.
 *           - expandAxis should be logically the same as taking the input function and copy/pasting it N times with different
 *             integers for each iteration (w/o the memory overhead of storing N exprtrees, of course)
 *             
 *              
 *             
 *           - how to "exit" a 'vary' loop is not clear
 *             - a node which converts a buffer to a function
 *               - what to do at boundaries?
 *             - early exit is a bit foggy too.
 *             - ergo while() is foggy. do we allocate infinity cells and then bail out early? hardly.
 *               - the way to do this is recursion.
 *               - if you need a long list, then you slowly construct a function of int which grows in scope,
 *                 then sample this via a buffer
 *                 - this doesn't conceptually correlate to the idea of "vary" nodes terribly well,
 *                   might want to provide a convenience wrapper
 *           - important: should *not* allocate cells until they are needed, e.g. if 
 *             we have nested vary()s to iterate upon multiple axes
 *             - e.g. we have an image which is a function of x,y. don't allocate an array of size X, only
 *               to discard it and allocate an array of size X*Y.
 *           - should be possible to build this function "in-place". almost like an anonymous function
 *           - what does this mean for vectors
 *           - we'd like to be able to treat a vector as a "thing". how do we do that with this scheme?
 *             - ability to package groups of data paths? i.e. x,y,z grouped together?
 *         - out of bounds behavior-- how is it handled?
 *           - can we just engineer things so it's impossible to go out of bounds?
 *             - i.e. functions are defined across all input, and iteration can't happen beyond the end of a buffer
 *         - how would we represent a matrix this way?
 *           - matrix : fn (ix,iy) -> v
 *         - it is then possible to compute axis coordinates the naive, expensive way, i.e. 
 *           with integer division and modulo.
 *           - well, duh. we are designing a programming language. can't stop the programmer from being
 *             totally stupid if they really want to. 
 *           - should just be intuitive to do it the smaert way 
 *           
 *   - what bucket-based evaluation gets us:
 *     - amortized cost of function call stack pushes
 *     - one-time evaluation of uniform inputs is possible
 *     - one-time type check
 *     - cache coherence? 
 *     - possibility to collapse identical inputs onto each other, to reduce re-evaluation.
 *     - parallelization becomes stupid easy, for:
 *       - CUDA/OpenCL implementation (beware of FP width!)
 *       - CPU vectorization
 *       - threading
 *   - MUTABLE (LARGE) OBJECTS:
 *     - example: a dynamically created expression tree
 *     - how about this:
 *       - buffers are updated if no one else depends on the previous state; copied if they do
 *     - do we need any other kind of mutable object?
 *     - need to be sure about order of operations-- with multithreading, we must be certain that an
 *       object-altering operation is complete in whatever thread it belongs to before the next
 *       operation which depends on it begins.
 *     - should we allow data to be added to expr nodes?
 *     - should we nodes of type "data"?
 *   - CONSTANCY SIMPLIFICATION:
 *     - know I am constant if both my inputs are constant
 *     - shared "token" which represents a cached value?
 *     - need to store a cache node everywhere a constant expr feeds to a non-constant expr.
 *       - otherwise that contributing branch will get re-evaluated
 *     - otherwise, cached value stored "beneath" const expr.
 *     - if all recipients of a const expr are also const, then caching is redundant.
 *     - perhaps this should all be rolled into the generic expr simplifier.
 *     - basically we are finding subtrees that have all-void leaves.
 *     - constancy needs to be understood across axes
 *       - i.e. a branch is constant within the context of an expand() evaluation
 *     - a procedure needs to be able to:
 *       - be used with multiple constancy "shapes"
 *   * BUFFER OWNERSHIP
 *     - when do we have to control buffer allocation?
 *       - when evaluating conditionals, we must split into True and False branch buffers
 *       - when evaluating expand()s, we must evaluate the expand length and create a 
 *         (flattened?) buffer for the new axes.
 *       - when evaluating a blocksize node
 *       - iterated function application-- when dest buffer is same as source.  
 *     - Where does the const expr short-circuit happen?
 *     - in Expr.eval(buffer,buffer) we have no knowledge of parents.
 *       - that's sort of not true anymore-- the Call node knows about its parents and could 
 *         allocate a buffer. The Expr node does not, but can receive the allocated buffer
 *         from Call. 
 *     - perhaps an Expr should return a newly allocated buffer, of size N or 1?
 *       - not great for buffer re-use; sharing.
 *       - however, we must ask expr to allocate output buffer at some point.
 *       - create nodes which alter the eval order and block size?
 *   - EDUCATION:
 *     - this language should be introduced as a game
 *     - user must match an output function using a limited selection/number of nodes
 *     - concepts and nodes are introduced slowly, one by one (like portal).
 *     - when a concept is introduced, it is the simplest possible example of a usage.
 *     - examples:
 *       - game to replace iterated loops with functional ones (reach a target exec speed?)
 *   
 *   - MAIN ELEGANCE ISSUES:
 *     - buffer ownership
 *       - when do we release a buffer (including effects of constant axes, blocks, sparseness, and recursion)
 *         > when we leave a function, we free the buffers.
 *         > this offers some control, and is a nice easy-to-understand concept: chunking.
 *     - determining evaluation order
 *       - including sparseness, multiple dependency, and recursion.
 *     - no control over evaluation order/buffer use. 
 *       - if the compiler decides to re-evaluate stupid things, we have no say over this.
 *       - do we have a way to tell what's being re-eval'd from userland?
 *       - this breaks reflection (?).
 *       - we may store things that we would rather re-compute on demand.
 *         and vice versa.
 *     - how do functions ask for/expect/specify input shape?
 *       - how strict/interpretive are we?
 *     - how do we avoid having to think about correlating data shape when we want it to mean
 *       iteration over distinct parts instead?
 *       - basically, I don't want to be forced to visualize "axes" and "dimensions" when all I want
 *         to do is iterate over some stuff.
 *     - no way to specify breadth or depth
 *       - do we evaluate the entire domain within a function, or do we descend through it in chunks?
 *       - breadth is speedy, depth is memory-light.
 *       - should be able to do it either way for different parts of the program.
 *       - should even be able to dynamically choose. i.e. argument to a function. 
 *     - basically, the evaluation system is much too black-boxy. 
 *     
 *     I FUCKING LOVE
 *       - how classes come naturally.
 *         - we just define a new object with const functions. the compiler can
 *           prove the objects aren't dynamic and will just dereference them cheaply.
 *         - you want virtual functions? then just make your functions varying.
 */

#ifndef DYNAMICEXPR_H_
#define DYNAMICEXPR_H_

#include <typeinfo>
#include <stdexcept>
#include <boost/intrusive_ptr.hpp>

namespace geom {
/*
class type_mismatch : public std::runtime_error {
public:
    type_mismatch(const std::string &msg):std::runtime_error(msg){}
    virtual ~type_mismatch() throw (){};
};

class buffer_mismatch : public std::runtime_error {
    buffer_mismatch(const std::string &msg):std::runtime_error(msg){}
    virtual ~type_mismatch() throw (){};
};

using boost::intrusive_ptr;

class DynamicExpr;

typedef boost::intrusive_ptr<DynamicExpr> DynExprPtr;

*/


/***************************
 * Dynamic expr class      *
 ***************************/

/*
class DynamicExpr {
    
    DynamicExpr():refct(0){}
    
    virtual void eval(DynTypePtr out, DynTypePtr in) = 0;
    virtual void eval(DynTypeArray out, DynTypeArray in) = 0;
    virtual vector<DynExprPtr> getInputs() = 0;
    virtual const type_info& outputType() = 0;
    virtual const type_info& inputType() = 0;
    
    
private:
    friend void intrusive_ptr_add_ref(DynamicExpr *pExpr);
    friend void intrusive_ptr_release(DynamicExpr *pExpr);
    
    size_t refct;
};
*/

/***************************
 * Reference counting      *
 ***************************/

/*
void intrusive_ptr_add_ref(DynamicExpr *pExpr){
    pExpr->refct++;
}

void intrusive_ptr_release(DynamicExpr *pExpr){
    if ((pExpr->refct -= 1) <= 0){
        delete pExpr;
    }
}
*/

}; //end namespace geom

#endif /* DYNAMICEXPR_H_ */
