/*
 * Expression.hpp
 *
 *  Created on: Feb 21, 2010
 *      Author: tbabb
 */


#ifndef EXPRESSION_HPP_
#define EXPRESSION_HPP_

#include <boost/shared_ptr.hpp>
#include <string>
#include <cmath>
#include <iostream>

#include <geomc/function/DynamicExpr.h>
#include <geomc/function/FunctionTypes.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/AffineTransform.h>

//macros: some shorthands for common verbose things
#define ExprPtr(I,O) boost::shared_ptr< Expr< I,O > >
#define ExprKindPtr(K,I,O) boost::shared_ptr< K< I,O > > //e.g. ptr(ConstExpr<double,int>) == ExprKindPtr(ConstExpr,double,int)
#define TEMPL1(A) template <typename A> 
#define TEMPL2(A,B) template <typename A, typename B>
#define TEMPL3(A,B,C) template <typename A, typename B, typename C>
#define TEMPL4(A,B,C,D) template <typename A, typename B, typename C, typename D>

namespace geom {

    using std::string;

    /*================================*
     * Wrapper Functions              *
     *================================*/

    // functions for use by binops
    TEMPL3(A,B,O) O add_o(A a, B b){ return a + b; }
    TEMPL3(A,B,O) O mul_o(A a, B b){ return a * b; }
    TEMPL3(A,B,O) O  lt_o(A a, B b){ return a < b; }
    TEMPL3(A,B,O) O  gt_o(A a, B b){ return a > b; }
    TEMPL3(A,B,O) O  eq_o(A a, B b){ return a == b; }
    TEMPL3(A,B,O) O leq_o(A a, B b){ return a <= b; }
    TEMPL3(A,B,O) O geq_o(A a, B b){ return a >= b; }
    TEMPL3(A,B,O) O min_o(A a, B b){ return a<b?a:b; }
    TEMPL3(A,B,O) O max_o(A a, B b){ return a>b?a:b; }
    
    //other functions
    TEMPL1(I) I null_fn(I i){ return i; }
    
    template <typename T, index_t N, index_t M>
    Vec<T,N> vec_reorder(Vec<T,M> sample, Vec<index_t,M> idxs){
        Vec<T,N> out;
        for (index_t i = 0; i < M; i++){
            out[idxs[i]] = sample[i];
        }
        return out;
    }
    
    template <typename T, index_t N>
    Vec<T,N+1> vec_collapse(Vec<T,N> input, T result){
        return Vec<T,N+1>(input, result);
    }

    /*================================*
     * Global Functions               *
     *================================*/

    //global function declarations
    TEMPL2(I,O)     ExprPtr(I,O) expr(O v);
    TEMPL1(I)       ExprPtr(I,I) nullx();
    TEMPL4(I,O,A,B) ExprPtr(I,O) binopx(ExprPtr(I,A) a, ExprPtr(I,B) b, O (*func)(A,B), string opname);
    TEMPL3(I,O,E)   ExprPtr(I,O) unopx(ExprPtr(I,E) a, O (*func)(E), string opname);
    TEMPL3(I,O,E)   ExprPtr(I,O) compose(ExprPtr(I,E) input_fn, ExprPtr(E,O) output_fn);

    TEMPL3(I,O,E)   ExprPtr(I,O) mixx(ExprPtr(I,E) k, ExprPtr(I,O) a, ExprPtr(I,O) b);
    TEMPL2(I,O)     ExprPtr(I,O) minx(ExprPtr(I,O) a, ExprPtr(I,O) b);
    TEMPL2(I,O)     ExprPtr(I,O) maxx(ExprPtr(I,O) a, ExprPtr(I,O) b);
    
    
    template <typename T, index_t N, typename O> 
    boost::shared_ptr< Expr<Vec<T,N>,O> > transformx(boost::shared_ptr< Expr<Vec<T,N>,O> > e, AffineTransform<T,N> xf);
    
    //like a "reorder", with the "order" varying over the input space
    template <typename T, index_t N, index_t M, typename O>
    boost::shared_ptr< Expr<Vec<T,M>,O> > slicex(boost::shared_ptr< Expr<Vec<T,N>,O> > sliced, boost::shared_ptr< Expr<Vec<T,M>,Vec<index_t,M> > > slicefn);
    
    //like a "reorder"
    template <typename T, index_t N, index_t M, typename O>
    boost::shared_ptr< Expr<Vec<T,M>,O> > slicex(boost::shared_ptr< Expr<Vec<T,N>,O> > sliced, Vec<index_t,M> sliceorder);
    
    //stack a bunch of exprs into the components of a Vec
    template <typename I, typename T, index_t N>
    boost::shared_ptr< Expr<I, Vec<T,N> > > stackx(boost::shared_ptr< Expr<I,T> > stack_fns[N]);
    
    //e.g., create a heightfield. Perlin2d -> Vec3d(x,y,perlin(x,y))
    //there is probably a concise mathematical term for this operation, but I don't know it.
    //Rn -> R => Rn -> R(n+1)
    template <typename T, index_t N>
    boost::shared_ptr< Expr<Vec<T,N>, Vec<T,N+1> > > hypersurfacex(boost::shared_ptr<Expr<Vec<T,N>,T> > surface);
    

    /*================================*
     * Expr Class                     *
     *================================*/

    template <typename I, typename O> class Expr { //: public DynamicExpr{
    public:

        typedef I in_type;
        typedef O out_type;

        Expr() {}
        virtual ~Expr() {}

        operator ExprPtr(I,O)(){ //casting
            return this->managedCopy();
        }
        
        /*
        void eval(DynTypeArray out, DynTypeArray in){
            //TODO: reduntant type check. we should just try to cast right away, and catch
            //if necessary.
            if (out.getElementType() != typeid(O) || in.getElementType() != typeid(I)){
                throw type_mismatch("buffer type mismatches input/output type of " + this->opname() + "\n" + \
                                        "function (I/O) : " + string(typeid(I).name()) + ", " + string(typeid(O).name()) + "\n" + \
                                        "buffer (I/O): " + string(typeid(I).name()) + ", " + string(typeid(O).name()) );
            } else if (in.size() != out.size() && in.size() != 1){
                 throw buffer_mismatch("buffer input/output sizes do not match");
            } else {
                O *outbuf = out.get<O>();
                if (in.size() == 1){
                    I &input = in.item<I>(0);
                    for (int i = 0; i < out.size(); i++){
                        outbuf[i] = eval(input);
                    }
                } else {
                    for (int i = 0; i < out.size(); i++){
                        
                    }
                }
            }
        }
        
        void eval(DynTypePtr out, DynTypePtr in){
            
        }*/

        ////// PURE VIRTUAL FUNCTIONS //////

        virtual O eval(I i) = 0;
        virtual void sout(std::ostream &s) = 0;
        virtual ExprPtr(I,O) managedCopy() = 0; //every individual subclass must override this!
        virtual string opname() = 0;
    };

    template <typename I, typename O>
    class ConstExpr:public Expr<I,O>{
    public:
        O v;

        ConstExpr():v(0){}
        ConstExpr(O value):v(value){};
        virtual ~ConstExpr(){};

        virtual void sout(std::ostream &s){
            s << v;
        }

        ConstExpr& operator+=(O d){
            v += d;
            return *this;
        }
        ConstExpr& operator*=(O d){
            v *= d;
            return *this;
        }

        virtual O eval(I x){
            return v;
        }

        virtual ExprPtr(I,O) managedCopy(){
            return ExprPtr(I,O)(new ConstExpr(v));
        }
        
        virtual string opname() { return "const"; }
    };
    
    //TODO: should be a singleton
    template <typename I>
    class NullExpr:public Expr<I,I>{
    public:

        NullExpr(){}
        virtual ~NullExpr(){};

        virtual void sout(std::ostream &s){
            s << "(null)";
        }

        virtual I eval(I x){
            return x;
        }

        virtual ExprPtr(I,I) managedCopy(){
            return ExprPtr(I,I)(new NullExpr());
        }
        
        virtual string opname() { return "null"; }
    };

    template <typename I, typename O, typename E>
    class MixExpr:public Expr<I,O>{
    public:
        ExprPtr(I,O) e1;
        ExprPtr(I,O) e2;
        ExprPtr(I,E) key;

        MixExpr(ExprPtr(I,E) key, ExprPtr(I,O) a, ExprPtr(I,O) b):e1(a),e2(b),key(key){}
        virtual ~MixExpr(){}

        virtual void sout(std::ostream &s){
            s << "mix(" << key << " , " << e1 << " : " << e2 <<")";
        }

        virtual O eval(I x){
            E s = key->eval(x);
            return (1-s)*e1->eval(x) + s*e2->eval(x);
        }

        virtual ExprPtr(I,O) managedCopy(){
            return ExprPtr(I,O)(new MixExpr(*this));
        }
        
        virtual string opname() { return "mix"; }
    };

    template <typename I, typename O, typename E>
    class ComposedExpr:public Expr<I,O>{
    public:
        ExprPtr(I,E) inner; //i.e. g(x) in f(g(x))
        ExprPtr(E,O) outer;

        ComposedExpr(ExprPtr(I,E) inner, ExprPtr(E,O) outer):inner(inner),outer(outer){}
        virtual ~ComposedExpr(){}

        virtual void sout(std::ostream &s){
            s << "(" << inner << " --> " << outer << ")";
        }

        virtual O eval(I x){
            return outer->eval(inner->eval(x));
        }

        virtual ExprPtr(I,O) managedCopy(){
            return ExprPtr(I,O)(new ComposedExpr<I,O,E>(*this));
        }
        
        virtual string opname() { return "compose"; }
    };

    template <typename I, typename O, typename A, typename B>
    class BinopExpr:public Expr<I,O>{
    public:
        ExprPtr(I,A) e1;
        ExprPtr(I,B) e2;
        O (*combine)(A, B);
        string name;

        BinopExpr(ExprPtr(I,A) a, ExprPtr(I,B) b, O (*combine)(A,B), string name):
            e1(a),e2(b),combine(combine),name(name){}
        virtual ~BinopExpr(){}

        virtual void sout(std::ostream &s){
            s << "(" << e1 << " " << (this->name) << " " << e2 << ")";
        }

        virtual O eval(I x){
            return (*combine)(e1->eval(x), e2->eval(x));
        }

        virtual ExprPtr(I,O) managedCopy(){
            return ExprPtr(I,O)(new BinopExpr<I,O,A,B>(*this));
        }
        
        virtual string opname() { return name; }
    };

    template <typename I, typename O, typename E>
    class UnaryOpExpr:public Expr<I,O>{
    public:
        ExprPtr(I,E) f;
        O (*func)(E);
        string name;

        UnaryOpExpr(ExprPtr(I,E) f, O (*func)(E), string name):
            f(f),func(func),name(name){}
        virtual ~UnaryOpExpr(){}

        virtual void sout(std::ostream &s){
            s << this->name << "(" << this->f << ")";
        }

        virtual O eval(I x){
            return (*func)(f->eval(x));
        }

        virtual ExprPtr(I,O) managedCopy(){
            return ExprPtr(I,O)(new UnaryOpExpr(*this));
        }
        virtual string opname() { return name; }
    };
    
    template <typename I, typename T, index_t N> 
    class StackedExpr : public Expr <I, Vec<T,N> > {
    protected:
        boost::shared_ptr< Expr<I,T> > fns[N];
    public:
        
        StackedExpr(boost::shared_ptr< Expr<I,T> > stackd_fns[N]){
            for (index_t i = 0; i < N; i++){
                fns[i] = stackd_fns[i];
            }
        }
        
        virtual ~StackedExpr(){}
        
        virtual void sout(std::ostream &s){
            s << "stack(";
            for (index_t i = 0; i < N; i++){
                s << fns[i];
                if (i < N-1){
                    s << ", ";
                }
            }
            s << ")";
        }
        
        virtual Vec<T,N> eval(I x){
            Vec<T,N> o;
            for (index_t i = 0; i < N; i++){
                o[i] = fns[i]->eval(x);
            }
            return o;
        }
        
        virtual boost::shared_ptr< Expr<I, Vec<T,N> > > managedCopy(){
            return boost::shared_ptr< Expr<I, Vec<T,N> > >(new StackedExpr(*this));
        }
        
        virtual string opname() { return "stack"; }
    };
    
    template <typename T, index_t N, typename O>
    class TransformedExpr : public Expr<Vec<T,N>, O>{
    public:
        boost::shared_ptr< Expr<Vec<T,N>,O> > e;
        AffineTransform<T,N> xf;
        
        TransformedExpr(boost::shared_ptr< Expr<Vec<T,N>,O> > e, AffineTransform<T,N> xf):
            e(e),xf(xf){}
        virtual ~TransformedExpr(){}
        
        virtual void sout(std::ostream &s){
            s << "transform(" << this->e << ")";
        }
        
        virtual O eval(Vec<T,N> x){
            return e->eval(x / xf);
        }
        
        virtual boost::shared_ptr< Expr<Vec<T,N>,O> > managedCopy(){
            return boost::shared_ptr< Expr<Vec<T,N>,O> >(new TransformedExpr(*this));
        }
        
        virtual string opname() { return "transform"; }
    };
    
    /*************************************
     * Wrapper Exprs                     *
     *************************************/
    
    template <typename T, index_t N>
    class PathExpr : public Expr< T, Vec<T,N> > {
    public:
        Path<T,N> path;
        
        PathExpr() {}
        PathExpr(const Path<T,N> &p):path(p) {}
        virtual ~PathExpr() {}
        
        virtual void sout(std::ostream &s) {
            s << "path(";
            
            for (size_t i = 0; i < path.knots.size(); i++){
                s << path.knots[i] << ", ";
            }
            
            s << ")";
        }
        
        virtual Vec<T,N> eval(T x) {
            return path.eval(x);
        }
        
        virtual boost::shared_ptr< Expr<T, Vec<T,N> > > managedCopy() {
            return boost::shared_ptr< Expr<T, Vec<T,N> > >(new PathExpr(*this));
        }
        
        virtual string opname() { return "path"; }
    };
    
    template <typename T, index_t N>
    class PerlinNoiseExpr : public Expr< Vec<T,N>, T > {
    public:
        PerlinNoise<T,N> perlin;
        
        PerlinNoiseExpr() {}
        PerlinNoiseExpr(const PerlinNoise<T,N> &p):perlin(p) {}
        virtual ~PerlinNoiseExpr() {}
        
        virtual void sout(std::ostream &s) {
            s << "perlin()";
        }
        
        virtual T eval(Vec<T,N> x) {
            return perlin.eval(x);
        }
        
        virtual boost::shared_ptr< Expr<Vec<T,N>,T > > managedCopy() {
            return boost::shared_ptr< Expr<Vec<T,N>,T > >(new PerlinNoiseExpr(*this));
        }
        
        virtual string opname() { return "perlin"; }
    };

    /*************************************
     * Template Function Bodies          *
     *************************************/

    template <typename I, typename O> ExprPtr(I,O) expr(O v){
        return ExprPtr(I,O)(new ConstExpr<I,O>(v));
    }
    
    template <typename I, typename O, typename A, typename B> ExprPtr(I,O) binopx(ExprPtr(I,A) a, ExprPtr(I,B) b, O (*func)(A,B), string opname){
        ExprKindPtr(ConstExpr,I,A) c1;
        ExprKindPtr(ConstExpr,I,B) c2;
        if ((c1 = boost::dynamic_pointer_cast<ConstExpr<I,A> >(a)) && (c2 = boost::dynamic_pointer_cast<ConstExpr<I,B> >(b))){
            return expr<I,O>((*func)(c1->v,c2->v));
        } else {
            return ExprPtr(I,O)(new BinopExpr<I,O,A,B>(a,b,func,opname));
        }
    }
    
    template <typename I, typename O, typename E> ExprPtr(I,O) unopx(ExprPtr(I,E) a, O (*func)(E), string opname){
        ExprKindPtr(ConstExpr,I,E) c;
        if ((c = boost::dynamic_pointer_cast<ConstExpr<I,E> >(a)) != NULL){
            return expr<I,O>((*func)(c->v));
        } else {
            return ExprPtr(I,O)(new UnaryOpExpr<I,O,E>(a,func,opname));
        }
    }
    
    TEMPL1(I) ExprPtr(I,I) nullx(){
        return ExprPtr(I,I)(new NullExpr<I>());
    }
    
    TEMPL3(I,O,E) ExprPtr(I,O) compose(ExprPtr(I,E) input_fn, ExprPtr(E,O) output_fn){
        return ExprPtr(I,O)(new ComposedExpr<I,O,E>(input_fn,output_fn));
    }
    
    TEMPL3(I,O,E) ExprPtr(I,O) mixx(ExprPtr(I,E) k, ExprPtr(I,O) a, ExprPtr(I,O) b){
        return MixExpr<I,O,E>(k,a,b);
    }
    
    TEMPL2(I,O) ExprPtr(I,O) minx(ExprPtr(I,O) a, ExprPtr(I,O) b){
        return binopx(a, b, &min_o<O,O,O>, "min");
    }
    
    TEMPL2(I,O) ExprPtr(I,O) maxx(ExprPtr(I,O) a, ExprPtr(I,O) b){
        return binopx(a, b, &max_o<O,O,O>, "max");
    }
    
    template <typename T, index_t N, typename O> 
    boost::shared_ptr< Expr<Vec<T,N>,O> > transformx(boost::shared_ptr< Expr<Vec<T,N>,O> > e, AffineTransform<T,N> xf){
        return TransformedExpr<T,N,O>(e,xf);
    }
    
    //sliced is an O-valued function over N dimensions.
    //return an O-valued function over M dimensions, by remapping the axes of the input.
    //slicefn is a function over M dimensions which describes which axes to select
    template <typename T, index_t N, index_t M, typename O>
    boost::shared_ptr< Expr<Vec<T,M>,O> > slicex(boost::shared_ptr< Expr<Vec<T,N>,O> > sliced, boost::shared_ptr< Expr<Vec<T,M>, Vec<index_t,M> > > slicefn){
        boost::shared_ptr< Expr< Vec<T,M>, Vec<T,M> > > input = nullx< Vec<T,M> >();
        boost::shared_ptr< Expr< Vec<T,M>, Vec<T,N> > > shuffler = binopx(input, slicefn, vec_reorder<T,N,M>, "vector_reorder");
        return compose(shuffler, sliced);
    }
    
    template <typename T, index_t N, index_t M, typename O>
    boost::shared_ptr< Expr<Vec<T,M>,O> > slicex(boost::shared_ptr< Expr<Vec<T,N>,O> > sliced, Vec<index_t,M> sliceorder){
        return slicex(sliced, expr< Vec<T,M>,Vec<index_t,M> >(sliceorder));
    }
    
    template <typename I, typename T, index_t N>
    boost::shared_ptr< Expr<I, Vec<T,N> > > stackx(boost::shared_ptr< Expr<I,T> > stack_fns[N]){
        return boost::shared_ptr< Expr<I, Vec<T,N> > >(new StackedExpr<I,T,N>(stack_fns));
    }

    template <typename T, index_t N>
    boost::shared_ptr< Expr<Vec<T,N>, Vec<T,N+1> > > hypersurfacex(boost::shared_ptr<Expr<Vec<T,N>,T> > surface){
        return binopx(nullx<Vec<T,N> >(), surface, &(vec_collapse<T,N>), "hypersurface");
    }

    /*************************************
     * Operator Function Bodies          *
     *************************************/

    // + operator
    TEMPL3(I,O,E) ExprPtr(I,O) operator+(ExprPtr(I,O) a, ExprPtr(I,E) b){
        return binopx(a,b,&add_o<O,E,O>, "+");
    }

    TEMPL3(I,O,E) ExprPtr(I,O) operator+(ExprPtr(I,O) x, E v){
        return binopx(x,expr<I,E>(v), &add_o<O,E,O>, "+");
    }

    TEMPL3(I,O,E) ExprPtr(I,O) operator+(E v, ExprPtr(I,O) x){
        return binopx(x,expr<I,E>(v), &add_o<O,E,O>, "+");
    }

    // * operator
    TEMPL3(I,O,E) ExprPtr(I,O) operator*(ExprPtr(I,O) x, E v){
        return binopx(x,expr<I,E>(v), &mul_o<O,E,O>, "*");
    }

    TEMPL3(I,O,E) ExprPtr(I,O) operator*(E v, ExprPtr(I,O) x){
        return binopx(expr<I,E>(v),x, &mul_o<E,O,O>, "*");
    }

    TEMPL2(I,O) ExprPtr(I,O) operator*(ExprPtr(I,O) a, ExprPtr(I,O) b){
        return binopx(a,b,&mul_o<O,O,O>,"*");
    }

    // << stream operator
    TEMPL2(I,O) std::ostream& operator<< (std::ostream &s, ExprPtr(I,O) f){
        f->sout(s);
        return s;
    }

    //QUESTION: Functions w. more than one input? currying?

    //TODO: shared_ptr<const LinearField> X = NewLinearField(X_AXIS); //EDIT: HOW TO CREATE A CONST?

    //TODO: <, >, <=, >=, ==, maxx, minx, absx
    //TODO: sinx, cosx, tanx, asinx, acosx, atanx, atan2x
    //TODO: DistSquaredField(Vector3d pt)
    //TODO: PerlinNoise<Itype,dimensions> : public Expression<Itype, double>
    //TODO: FractalNoiseField()
    //TODO: RBFInterpField() ?
    //TODO: SplineNd (x --> Vector3d, e.g.)
    //TODO: Vector<Expression<double,double>, N) ?
    //TODO: Image: Expression< Vector<int,2>,Vector<double,3 or 4> >
    //TODO: InterpDiscrete : (vector int,N->value) -> (vector double,N->value)
    //TODO: FluidVolume: Expression<double,Vector3d>
    //TODO: OpticalFlow: (vector int,N->color) -> (somehow: int,int,int,...,double -> color) ??
    //TODO: Derivative: vector->double : vector->vector
    //TODO: Cache: define a region over which to store lazily evaluated values on a discrete space
            //fall back to a hashmap outside of region?


    //TODO: VoxelCacheField(field, sampleBox, gridsize)
    //TODO: BlurField() //should reference an internal CachedField to avoid extreme redundancy

} //end namespace geom

#endif /* EXPRESSION_HPP_ */
