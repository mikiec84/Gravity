//
//  func.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 24/10/16.
//
//
#include <cmath>
#include <gravity/func.h>
//
//
using namespace std;
namespace gravity{
  
    pair<func<double>,func<double>> get_mag_ang(const func<Cpx>& f){
        func<double> f_mag = get_mag(f._cst.get()), f_ang = get_ang(f._cst.get());
        for (auto &qt: *f._qterms) {
            auto coef = qt.second._coef;
            func<double> p1_mag, p2_mag, p1_ang, p2_ang;
            func<double> c_mag, c_ang;
            if(!qt.second._p->first->_polar){
                auto p1_r = get_real(qt.second._p->first.get());
                auto p1_i = get_imag(qt.second._p->first.get());
                p1_mag = sqrt(p1_r*p1_r + p1_i*p1_i);
                p1_ang = atan2(p1_r,p1_i);// TODO fix atan2
            }
            else {
                p1_mag = get_mag(qt.second._p->first.get());
                p1_ang = get_ang(qt.second._p->first.get());
            }
            if(!qt.second._p->second->_polar){
                auto p2_r = get_real(qt.second._p->second.get());
                auto p2_i = get_imag(qt.second._p->second.get());
                p2_mag = sqrt(p2_r*p2_r + p2_i*p2_i);
                p2_ang = atan2(p2_r,p2_i);// TODO fix atan2
            }
            else {
                p2_mag = get_mag(qt.second._p->second.get());
                p2_ang = get_ang(qt.second._p->second.get());
            }
            if(!coef->_polar){
                auto c_r = get_real(coef.get());
                auto c_i = get_imag(coef.get());
                c_mag = sqrt(c_r*c_r + c_i*c_i);
                c_ang = atan2(c_r,c_i);// TODO fix atan2
            }
            else {
                c_mag = get_mag(coef.get());
                c_ang = get_ang(coef.get());
            }
            if(coef->_is_transposed){
                c_mag.transpose();
                c_ang.transpose();
                p1_mag._is_vector = true;
                p1_ang._is_vector = true;
                p2_mag._is_vector = true;
                p2_ang._is_vector = true;
            }
            if(qt.second._sign){
                f_mag += c_mag*p1_mag*p2_mag;
                f_ang += c_ang + p1_ang + p2_ang;
            }
            else {
                f_mag -= c_mag*p1_mag*p2_mag;
                f_ang -= c_ang + p1_ang + p2_ang;
            }
        }
        for (auto &lt: *f._lterms) {
            auto coef = lt.second._coef;
            func<double> p1_mag, p1_ang;
            func<double> c_mag, c_ang;
            if(!lt.second._p->_polar){
                auto p1_r = get_real(lt.second._p.get());
                auto p1_i = get_imag(lt.second._p.get());
                p1_mag = sqrt(p1_r*p1_r + p1_i*p1_i);
                p1_ang = atan2(p1_r,p1_i);// TODO fix atan2
            }
            else {
                p1_mag = get_mag(lt.second._p.get());
                p1_ang = get_ang(lt.second._p.get());
            }
            if(!coef->_polar){
                auto c_r = get_real(coef.get());
                auto c_i = get_imag(coef.get());
                c_mag = sqrt(c_r*c_r + c_i*c_i);
                c_ang = atan2(c_r,c_i);// TODO fix atan2
            }
            else {
                c_mag = get_mag(coef.get());
                c_ang = get_ang(coef.get());
            }
            if(coef->_is_transposed){
                c_mag.transpose();
                c_ang.transpose();
                p1_mag._is_vector = true;
                p1_ang._is_vector = true;
            }
            if(lt.second._sign){
                f_mag += c_mag * p1_mag;
                f_ang += c_ang + p1_ang;
            }
            else {
                f_mag -= c_mag * p1_mag;
                f_ang -= c_ang + p1_ang;
            }
        }
        return {f_mag,f_ang};
    }
    
    /** WARNING, only call if the variables appearing in the function are complex or double */
    pair<func<double>,func<double>> get_real_imag(const func<Cpx>& f){
        func<double> c_r, c_i, c_mag, c_ang;
        func<double> p1_r, p2_r, p1_i, p2_i;
        func<double> p1_mag, p2_mag, p1_ang, p2_ang;
        func<double> f_r = get_real(f._cst.get()), f_i = get_imag(f._cst.get());
        for (auto &pt: *f._pterms) {
            func<> fp_r = 1.;
            func<> fp_i = 0.;
            auto coef = pt.second._coef;
            if(coef->_polar){
                c_mag = get_mag(coef.get());
                c_ang = get_ang(coef.get());
                c_r = c_mag*cos(c_ang);
                c_i = c_mag*sin(c_ang);
            }
            else {
                c_r = get_real(coef.get());
                c_i = get_imag(coef.get());
            }
            if(coef->_is_transposed){
                c_r.transpose();
                c_i.transpose();
                p1_r._is_vector = true;
                p1_i._is_vector = true;
                p2_r._is_vector = true;
                p2_i._is_vector = true;
            }
            auto it = pt.second._l->begin();
            while(it!=pt.second._l->end()){
                auto vv = it->first;
                auto expo = it->second;
                if(vv->_polar){
                    p1_mag = get_mag(vv.get());
                    p1_ang = get_ang(vv.get());
                    p1_r = p1_mag*cos(p1_ang);
                    p1_i = p1_mag*sin(p1_ang);
                }
                else {
                    p1_r = get_real(vv.get());
                    p1_i = get_imag(vv.get());
                }
                if(coef->_polar || vv->_polar){
                    throw invalid_argument("unsupported yet");
                }
                else {
                    for (auto i = 0; i<expo; i++) {
                        auto fp_r_new = fp_r*p1_r - fp_i*p1_i;
                        auto fp_i_new = fp_r*p1_i + fp_i*p1_r;
                        fp_r = fp_r_new;
                        fp_i = fp_i_new;
                    }
                }
                it++;
            }
            auto fp_r_new = c_r*fp_r - c_i*fp_i;
            auto fp_i_new = c_r*fp_i + c_i*fp_r;
            fp_r = fp_r_new;
            fp_i = fp_i_new;
            if(!pt.second._sign){
                fp_r.reverse_sign();
                fp_i.reverse_sign();
            }
            f_r += fp_r;
            f_i += fp_i;
        }
        for (auto &qt: *f._qterms) {
            auto coef = qt.second._coef;
            func<double> p1_r, p2_r, p1_i, p2_i;
            func<double> p1_mag, p2_mag, p1_ang, p2_ang;
            func<double> c_r, c_i, c_mag, c_ang;
            func<> fq_r, fq_i;
            if(qt.second._p->first->_polar){
                p1_mag = get_mag(qt.second._p->first.get());
                p1_ang = get_ang(qt.second._p->first.get());
                p1_r = p1_mag*cos(p1_ang);
                p1_i = p1_mag*sin(p1_ang);
            }
            else {
                p1_r = get_real(qt.second._p->first.get());
                p1_i = get_imag(qt.second._p->first.get());
            }
            if(qt.second._p->second->_polar){
                p2_mag = get_mag(qt.second._p->second.get());
                p2_ang = get_ang(qt.second._p->second.get());
                p2_r = p2_mag*cos(p2_ang);
                p2_i = p2_mag*sin(p2_ang);
            }
            else {
                p2_r = get_real(qt.second._p->second.get());
                p2_i = get_imag(qt.second._p->second.get());
            }
            if(qt.second._p->first->_is_transposed){
                p1_r.transpose();
                p1_i.transpose();
                p2_r._is_vector = true;
                p2_i._is_vector = true;
            }
            if(coef->_polar){
                c_mag = get_mag(coef.get());
                c_ang = get_ang(coef.get());
                c_r = c_mag*cos(c_ang);
                c_i = c_mag*sin(c_ang);
            }
            else {
                c_r = get_real(coef.get());
                c_i = get_imag(coef.get());
            }
            if(coef->_is_transposed){
                c_r.transpose();
                c_i.transpose();
                p1_r._is_vector = true;
                p1_i._is_vector = true;
                p2_r._is_vector = true;
                p2_i._is_vector = true;
            }
            if(coef->_polar && qt.second._p->first->_polar && qt.second._p->second->_polar){
                fq_r = c_mag*p1_mag*p2_mag*cos(c_ang+p1_ang+p2_ang);
                fq_i = c_mag*p1_mag*p2_mag*sin(c_ang+p1_ang+p2_ang);
            }
            else if(!coef->_polar && qt.second._p->first->_polar && qt.second._p->second->_polar){
                fq_r = c_r*p1_mag*p2_mag*cos(p1_ang+p2_ang) - c_i*p1_mag*p2_mag*sin(p1_ang+p2_ang);
                fq_i = c_r*p1_mag*p2_mag*sin(p1_ang+p2_ang) + c_i*p1_mag*p2_mag*cos(p1_ang+p2_ang);
            }
            else {
                fq_r = c_r*(p1_r*p2_r - p1_i*p2_i) - c_i*(p1_r*p2_i + p1_i*p2_r);
                fq_i = c_r*(p1_r*p2_i + p1_i*p2_r) + c_i*(p1_r*p2_r - p1_i*p2_i);
            }
            if(!qt.second._sign){
                fq_r.reverse_sign();
                fq_i.reverse_sign();
            }
            f_r += fq_r;
            f_i += fq_i;
        }
        for (auto &lt: *f._lterms) {
            auto coef = lt.second._coef;
            func<double> p1_r, p1_i, p1_mag, p1_ang;
            func<double> c_r, c_i, c_mag, c_ang;
            func<> fl_r, fl_i;
            if(lt.second._p->_polar){
                p1_mag = get_mag(lt.second._p.get());
                p1_ang = get_ang(lt.second._p.get());
                p1_r = p1_mag*cos(p1_ang);
                p1_i = p1_mag*sin(p1_ang);
            }
            else {
                p1_r = get_real(lt.second._p.get());
                p1_i = get_imag(lt.second._p.get());
            }
            if(coef->_polar){
                c_mag = get_mag(coef.get());
                c_ang = get_ang(coef.get());
                c_r = c_mag*cos(c_ang);
                c_i = c_mag*sin(c_ang);
            }
            else {
                c_r = get_real(coef.get());
                c_i = get_imag(coef.get());
            }
            if(coef->_is_transposed){
                c_r.transpose();
                c_i.transpose();
                p1_r._is_vector = true;
                p1_i._is_vector = true;
            }
            if(coef->_polar && lt.second._p->_polar){
                fl_r = c_mag*p1_mag*cos(c_ang+p1_ang);
                fl_i = c_mag*p1_mag*sin(c_ang+p1_ang);
            }
            else {
                fl_r = c_r*p1_r - c_i*p1_i;
                fl_i = c_r*p1_i + c_i*p1_r;
            }
            if(!lt.second._sign){
                fl_r.reverse_sign();
                fl_i.reverse_sign();
            }
            f_r += fl_r;
            f_i += fl_i;
        }
        return {f_r,f_i};
    }
//
//    bool is_indexed(const constant_* c){
//        if (c->is_var() || c->is_param()) {
//            return (((param_*)c)->is_indexed());
//        }
//        return true;
//    }
//
//    size_t get_id(const constant_* c){
//        if (c->is_var() || c->is_param()) {
//            return ((param_*)c)->get_vec_id();
//        }
//        return 0;
//    }
//
//    size_t get_id_inst(const constant_* c, unsigned inst){
//        if (c->is_var() || c->is_param()) {
//            return ((param_*)c)->get_id_inst(inst);
//        }
//        return 0;
//    }
//
////    void func_::eval_vector() {
////        //        try {
////        ////            if (_val->size()<_dim[0]) {
////        ////                _val->resize(_dim[0]);
////        ////            }
////        //            if (is_constant() && !_expr && _params->size()==1) { // This is a parameter
////        //                auto p_c = (*_params->begin()).second.first;
////        //                for (auto i = 0; i < p_c->get_dim(); i++) {
////        //                    set_val(i,t_eval(p_c.get(),i));
////        //                }
////        //                _evaluated = true;
////        //                return;
////        //            }
////        //        }
////        //        catch (const std::out_of_range& e) {
////        //            cout << "Out of Range before eval.";
////        //        }
////        //        if(_expr && _expr->is_bexpr()){
////        //            auto be = (bexpr*)_expr.get();
////        //            if (be->_otype==product_ && _lterms->empty() && _pterms->empty() && _qterms->empty()) {//Pure product
////        //                if (be->_lson->is_matrix() && be->_rson->_is_vector) {//matrix*vect
////        //                    for (size_t i = 0; i < _dim[0]; i++) {
////        //                        double res = 0;
////        //                        for (size_t j = 0; j < be->_lson->_dim[1]; j++) {
////        //                            res += be->_lson->get_val(i,j) * be->_rson->get_val(j);
////        //                        }
////        //                        set_val(i,be->_coef*res);
////        //                    }
////        //                    return;
////        //                }
////        //            }
////        //        }
////        if (_is_transposed) {
////            for (auto inst = 0; inst < _dim[1]; inst++) {
////                set_val(inst,eval(inst));
////            }
////        }
////        else {
////            for (auto inst = 0; inst < _dim[0]; inst++) {
////                set_val(inst,eval(inst));
////            }
////        }
////        //        }
////        //        catch (const std::out_of_range& e) {
////        //            cout << "Out of Range error at inst = " << inst << endl;
////        //        }
////
////        //        }
////    }
//
//
////    void func_::eval_matrix() {
////        //        if (is_constant() && !_expr && _params->size()==1) { // This is a parameter
////        //            auto p_c = (*_params->begin()).second.first;
////        //            for (size_t i = 0; i < _dim[0]; i++) {
////        //                for (size_t j = 0; j < _dim[1]; j++) {
////        //                    set_val(i,j,t_eval(p_c.get(), i, j));
////        //                }
////        //            }
////        //            _evaluated = true;
////        //            return;
////        //        }
////        //        double res = 0;
////        if(!_expr){
////            for (size_t i = 0; i < _dim[0]; i++) {
////                for (size_t j = 0; j < _dim[1]; j++) {
////                    set_val(i,j,eval(i,j));
////                }
////            }
////        }
////        if(is_constant()){
////            _evaluated = true;
////        }
////        if(_expr){
////            auto be = (bexpr*)_expr.get();
////            if (be->_otype==product_ && _lterms->empty() && _pterms->empty() && _qterms->empty()) {//Pure product
////                if (be->_lson->is_matrix() && be->_rson->is_matrix()) {
////                    //matrix product
////                    if (!_is_hessian) {
////                        for (size_t i = 0; i < _dim[0]; i++) {
////                            for (size_t j = 0; j < _dim[1]; j++) {
////                                auto res = 0.;
////                                for (size_t col = 0; col<be->_lson->_dim[1]; col++) {
////                                    res += be->_lson->get_val(i,col) * be->_rson->get_val(col,j);
////                                }
////                                set_val(i,j, be->_coef*res);
////                            }
////                        }
////                    }
////                    else {
////                        for (size_t i = 0; i < _dim[0]; i++) {
////                            for (size_t j = i; j < _dim[1]; j++) {
////                                auto res = 0.;
////                                for (size_t col = 0; col<be->_lson->_dim[1]; col++) {
////                                    res += be->_lson->get_val(i,col) * be->_rson->get_val(col,j);
////                                }
////                                set_val(i,j, be->_coef*res);
////                            }
////                        }
////                    }
////                    return;
////                }
////                if (be->_lson->is_matrix() && !be->_rson->is_matrix() && be->_rson->_is_transposed) {//matrix * transposed vect
////                    for (size_t i = 0; i < _dim[0]; i++) {
////                        for (size_t j = 0; j < _dim[1]; j++) {
////                            set_val(i,j, be->_coef*be->_lson->get_val(i,j) * be->_rson->get_val(j));
////                        }
////                    }
////                    return;
////                }
////                if (!be->_lson->is_matrix() && !be->_lson->_is_transposed && be->_rson->is_matrix() ) {//vect * matrix
////                    for (size_t i = 0; i < _dim[0]; i++) {
////                        for (size_t j = 0; j < _dim[1]; j++) {
////                            set_val(i,j, be->_coef*(be->_lson->get_val(i) * be->_rson->get_val(i,j)));
////                        }
////                    }
////                    return;
////                }
////                //                if (be->_lson->is_matrix() && be->_rson->_is_vector) {//matrix*vect
////                //                    for (size_t i = 0; i < _dim[0]; i++) {
////                //                        for (size_t j = 0; j < _dim[1]; j++) {
////                //                            set_val(i,j, be->_coef*(be->_lson->get_val(i,j) * be->_rson->get_val(i)));
////                //                        }
////                //                    }
////                //                    return;
////                //                }
////            }
////            if (!_is_hessian) {
////                for (size_t i = 0; i < _dim[0]; i++) {
////                    for (size_t j = 0; j < _dim[1]; j++) {
////                        auto res = 0.;
////                        for (size_t col = 0; col<be->_lson->_dim[1]; col++) {
////                            res += be->_lson->get_val(i,col) * be->_rson->get_val(col,j);
////                        }
////                        set_val(i,j,res);
////                    }
////                }
////            }
////            else {
////                for (size_t i = 0; i < _dim[0]; i++) {
////                    for (size_t j = i; j < _dim[1]; j++) {
////                        auto res = 0.;
////                        for (size_t col = 0; col<be->_lson->_dim[1]; col++) {
////                            res += be->_lson->get_val(i,col) * be->_rson->get_val(col,j);
////                        }
////                        set_val(i,j,res);
////                    }
////                }
////            }
////        }
////    }
//
////    double t_eval(const constant_* c, size_t i){
////        if (!c) {
////            throw invalid_argument("Cannot evaluate nullptr!");
////        }
////        switch (c->get_type()) {
////            case binary_c: {
////                return ((constant<bool>*)(c))->eval();
////                break;
////            }
////            case short_c: {
////                return ((constant<short>*)(c))->eval();
////                break;
////            }
////            case integer_c: {
////                return ((constant<int>*)(c))->eval();
////                break;
////            }
////            case float_c: {
////                return ((constant<float>*)(c))->eval();
////                break;
////            }
////            case double_c: {
////                return ((constant<double>*)(c))->eval();
////                break;
////            }
////            case long_c: {
////                return (double)((constant<long double>*)(c))->eval();
////                break;
////            }
////            case par_c:{
////                auto p_c2 = (param_*)(c);
////                switch (p_c2->get_intype()) {
////                    case binary_:
////                        return ((param<bool>*)p_c2)->eval(i);
////                        break;
////                    case short_:
////                        return ((param<short>*)p_c2)->eval(i);
////                        break;
////                    case integer_:
////                        return ((param<int>*)p_c2)->eval(i);
////                        break;
////                    case float_:
////                        return ((param<float>*)p_c2)->eval(i);
////                        break;
////                    case double_:
////                        return ((param<double>*)p_c2)->eval(i);
////                        break;
////                    case long_:
////                        return (double)((param<long double>*)p_c2)->eval(i);
////                        break;
////                    default:
////                        break;
////                }
////                break;
////            }
////            case var_c:{
////                auto p_c2 = (param_*)(c);
////                switch (p_c2->get_intype()) {
////                    case binary_:
////                        return ((var<bool>*)p_c2)->eval(i);
////                        break;
////                    case short_:
////                        return ((var<short>*)p_c2)->eval(i);
////                        break;
////                    case integer_:
////                        return ((var<int>*)p_c2)->eval(i);
////                        break;
////                    case float_:
////                        return ((var<float>*)p_c2)->eval(i);
////                        break;
////                    case double_:
////                        return ((var<double>*)p_c2)->eval(i);
////                        break;
////                    case long_:
////                        return (double)((var<long double>*)p_c2)->eval(i);
////                        break;
////                    default:
////                        break;
////                }
////                break;
////            }
////            case uexp_c: {
////                return ((uexpr*)c)->eval(i);
////                break;
////            }
////            case bexp_c: {
////                return ((bexpr*)c)->eval(i);
////                break;
////            }
////            case func_c: {
////                auto f = (func_*)c;
////                if (f->is_constant() && f->_evaluated) {
////                    if (f->is_number()) {
////                        return f->_val->at(0);
////                    }
////                    return f->get_val(i);
////                }
////                return ((func_*)c)->eval(i);
////                break;
////            }
////            default:
////                break;
////        }
////        return 0;
////    }
//
////    double t_eval(const constant_* c, size_t i, size_t j){
////        if (!c) {
////            throw invalid_argument("Cannot evaluate nullptr!");
////        }
////        //        if (!c->is_matrix()) {
////        //            return eval(c, j);
////        //        }
////        switch (c->get_type()) {
////            case binary_c: {
////                return ((constant<bool>*)(c))->eval();
////                break;
////            }
////            case short_c: {
////                return ((constant<short>*)(c))->eval();
////                break;
////            }
////            case integer_c: {
////                return ((constant<int>*)(c))->eval();
////                break;
////            }
////            case float_c: {
////                return ((constant<float>*)(c))->eval();
////                break;
////            }
////            case double_c: {
////                return ((constant<double>*)(c))->eval();
////                break;
////            }
////            case long_c: {
////                return (double)((constant<long double>*)(c))->eval();
////                break;
////            }
////            case par_c:{
////                auto p_c2 = (param_*)(c);
////                switch (p_c2->get_intype()) {
////                    case binary_:
////                        return ((param<bool>*)p_c2)->eval(i,j);
////                        break;
////                    case short_:
////                        return ((param<short>*)p_c2)->eval(i,j);
////                        break;
////                    case integer_:
////                        return ((param<int>*)p_c2)->eval(i,j);
////                        break;
////                    case float_:
////                        return ((param<float>*)p_c2)->eval(i,j);
////                        break;
////                    case double_:
////                        return ((param<double>*)p_c2)->eval(i,j);
////                        break;
////                    case long_:
////                        return (double)((param<long double>*)p_c2)->eval(i,j);
////                        break;
////                    default:
////                        break;
////                }
////                break;
////            }
////            case var_c:{
////                auto p_c2 = (param_*)(c);
////                switch (p_c2->get_intype()) {
////                    case binary_:
////                        return ((var<bool>*)p_c2)->eval(i,j);
////                        break;
////                    case short_:
////                        return ((var<short>*)p_c2)->eval(i,j);
////                        break;
////                    case integer_:
////                        return ((var<int>*)p_c2)->eval(i,j);
////                        break;
////                    case float_:
////                        return ((var<float>*)p_c2)->eval(i,j);
////                        break;
////                    case double_:
////                        return ((var<double>*)p_c2)->eval(i,j);
////                        break;
////                    case long_:
////                        return (double)((var<long double>*)p_c2)->eval(i,j);
////                        break;
////                    default:
////                        break;
////                }
////                break;
////            }
////            case uexp_c: {
////                return ((uexpr*)c)->eval(i,j);
////                break;
////            }
////            case bexp_c: {
////                return ((bexpr*)c)->eval(i,j);
////                break;
////            }
////            case func_c: {
////                auto f = (func_*)c;
////                if (f->is_constant() && f->_evaluated) {
////                    if (f->is_number()) {
////                        return f->_val->at(0);
////                    }
////                    if (f->is_matrix()) {
////                        return f->get_val(i, j);
////                    }
////                    return f->get_val(j);
////                }
////                return ((func_*)c)->eval(i,j);
////                break;
////            }
////            default:
////                break;
////        }
////        return 0;
////    }
//
//
//    func_ get_derivative(constant_* c, const param_ &v){
//        if (!c) {
//            throw invalid_argument("Cannot evaluate nullptr!");
//        }
//        if (c->is_number() || c->is_param()) {
//            return func_();
//        }
//        if (c->is_var()) {
//            if ((*(param_*)c)==v) {
//                func_() += 1;
//            }
//            return func_();
//        }
//
//        if(c->is_uexpr()){
//            return ((uexpr*)c)->get_derivative(v);
//        }
//        if(c->is_bexpr()){
//            return ((bexpr*)c)->get_derivative(v);
//        }
//        if(c->is_function()){
//            return ((func_*)c)->get_derivative(v);
//        }
//        return func_();
//    }
//
//
//    Sign constant_::get_all_sign() const{
//        switch (_type) {
//            case binary_c: {
//                return ((constant<bool>*)this)->get_sign();
//                break;
//            }
//            case short_c: {
//                return ((constant<short>*)this)->get_sign();
//                break;
//            }
//            case integer_c: {
//                return ((constant<int>*)this)->get_sign();
//                break;
//            }
//            case float_c: {
//                return ((constant<float>*)this)->get_sign();
//                break;
//            }
//            case double_c: {
//                return ((constant<double>*)this)->get_sign();
//                break;
//            }
//            case long_c: {
//                return ((constant<long double>*)this)->get_sign();
//                break;
//            }
//            case par_c:{
//                return ((param_*)this)->get_all_sign();
//                break;
//            }
//            case uexp_c: {
//                return ((uexpr*)this)->get_all_sign(); // TO UPDATE
//                break;
//            }
//            case bexp_c: {
//                return ((bexpr*)this)->get_all_sign(); // TO UPDATE
//                break;
//            }
//            case var_c:{
//                return ((param_*)this)->get_all_sign();
//                break;
//            }
//            case func_c: {
//                return ((func_*)this)->get_all_sign();
//                break;
//            }
//            default:
//                break;
//        }
//        return unknown_;
//    }
//
//
//
//    Sign constant_::get_sign(size_t idx) const{
//        switch (_type) {
//            case binary_c: {
//                return ((constant<bool>*)this)->get_sign();
//                break;
//            }
//            case short_c: {
//                return ((constant<short>*)this)->get_sign();
//                break;
//            }
//            case integer_c: {
//                return ((constant<int>*)this)->get_sign();
//                break;
//            }
//            case float_c: {
//                return ((constant<float>*)this)->get_sign();
//                break;
//            }
//            case double_c: {
//                return ((constant<double>*)this)->get_sign();
//                break;
//            }
//            case long_c: {
//                return ((constant<long double>*)this)->get_sign();
//                break;
//            }
//            case par_c:{
//                return ((param_*)this)->get_sign(idx);
//                break;
//            }
//            case uexp_c: {
//                return ((uexpr*)this)->get_sign(idx);
//                break;
//            }
//            case bexp_c: {
//                return ((bexpr*)this)->get_sign(idx);
//                break;
//            }
//            case var_c:{
//                return ((param_*)this)->get_sign(idx);
//                break;
//            }
//            case func_c: {
//                return ((func_*)this)->get_sign(idx);
//                break;
//            }
//            default:
//                break;
//        }
//        return unknown_;
//    }
//
//    Sign param_::get_sign(size_t idx) const{
//        switch (_intype) {
//            case binary_:
//                if (is_param()) {
//                    return ((param<bool>*)this)->get_sign(idx);
//                }
//                return ((var<bool>*)this)->get_sign(idx);
//                break;
//            case short_:
//                if (is_param()) {
//                    return ((param<short>*)this)->get_sign(idx);
//                }
//                return ((var<short>*)this)->get_sign(idx);
//                break;
//            case integer_:
//                if (is_param()) {
//                    return ((param<int>*)this)->get_sign(idx);
//                }
//                return ((var<int>*)this)->get_sign(idx);
//                break;
//            case float_:
//                if (is_param()) {
//                    return ((param<float>*)this)->get_sign(idx);
//                }
//                return ((var<float>*)this)->get_sign(idx);
//                break;
//            case double_:
//                if (is_param()) {
//                    return ((param<double>*)this)->get_sign(idx);
//                }
//                return ((var<double>*)this)->get_sign(idx);
//                break;
//            case long_:
//                if (is_param()) {
//                    return ((param<long double>*)this)->get_sign(idx);
//                }
//                return ((var<long double>*)this)->get_sign(idx);
//                break;
//            default:
//                throw invalid_argument("Unsupported type");
//                break;
//        }
//    }
//
//    Sign param_::get_all_sign() const{
//        switch (_intype) {
//            case binary_:
//                return pos_;
//                break;
//            case short_:
//                if (is_param()) {
//                    return ((param<short>*)this)->get_all_sign<short>();
//                }
//                return ((var<short>*)this)->get_all_sign<short>();
//                break;
//            case integer_:
//                if (is_param()) {
//                    return ((param<int>*)this)->get_all_sign<int>();
//                }
//                return ((var<int>*)this)->get_all_sign<int>();
//                break;
//            case float_:
//                if (is_param()) {
//                    return ((param<float>*)this)->get_all_sign<float>();
//                }
//                return ((var<float>*)this)->get_all_sign<float>();
//                break;
//            case double_:
//                if (is_param()) {
//                    return ((param<double>*)this)->get_all_sign<double>();
//                }
//                return ((var<double>*)this)->get_all_sign<double>();
//                break;
//            case long_:
//                if (is_param()) {
//                    return ((param<long double>*)this)->get_all_sign<long double>();
//                }
//                return ((var<long double>*)this)->get_all_sign<long double>();
//                break;
//            case complex_:
//                if (is_param()) {
//                    return ((param<Cpx>*)this)->get_all_sign_cpx();
//                }
//                return ((var<Cpx>*)this)->get_all_sign_cpx();
//                break;
//            default:
//                return unknown_;
//                //                throw invalid_argument("Unsupported type");
//                break;
//        }
//
//    }
//
//
//
//
//
//
//
//    bool func_::operator==(const func_& f) const{
//        if (_ftype!=f._ftype || _all_sign!=f._all_sign || _vars->size()!=f._vars->size() || _pterms->size() != f._pterms->size()|| _qterms->size() != f._qterms->size() || _lterms->size() != f._lterms->size()) {
//            return false;
//        }
//        if (this->to_str()!=f.to_str()) {
//            return false;
//        }
//        return true;
//    }
//
//    func_& func_::operator=(const func_& f){
//        //        _evaluated = f._evaluated;
//        _to_str = f._to_str;
//
//        delete _all_range;
//        if (_vars) {
//            _vars->clear();
//        }
//        if (_params){
//            _params->clear();
//        }
//        _dfdx->clear();
//        delete _range;
//        delete _convexity;
//        delete _sign;
//        delete _lterms;
//        delete _qterms;
//        delete _pterms;
//        delete _DAG;
//        delete _queue;
//        delete _cst;
//
//        set_type(func_c);
//        _lterms = new map<string, lterm>();
//        _qterms = new map<string, qterm>();
//        _pterms = new map<string, pterm>();
//        _DAG = new map<string, expr*>();
//        _queue = new deque<shared_ptr<expr>>();
//        _ftype = f._ftype;
//        _return_type = f._return_type;
//        _all_convexity = f._all_convexity;
//        _all_sign = f._all_sign;
//        _is_transposed = f._is_transposed;
//        _is_vector = f._is_vector;
//        //        _is_matrix = f._is_matrix;
//        _is_constraint = f._is_constraint;
//        _is_hessian = f._is_hessian;
//        _dim[0] = f._dim[0];
//        _dim[1] = f._dim[1];
//        _embedded = f._embedded;
//        _return_type = f._return_type;
//        _cst = copy(*f._cst);
//        if (_cst->is_function()) {
//            embed(*(func_*)_cst);
//        }
//
//        _all_range = new pair<double,double>(f._all_range->first, f._all_range->second);
//        _sign = nullptr;
//        _convexity = nullptr;
//        _range = nullptr;
//
//        if(f._sign){
//            _sign = new vector<Sign>();
//            *_sign = *f._sign;
//        }
//        if (f._convexity) {
//            _convexity = new vector<Convexity>();
//            *_convexity = *f._convexity;
//        }
//        if (f._range) {
//            _range= new vector<pair<double,double>>();
//            *_range = *f._range;// Make sure this creates new pointers inside each pair, otherwise use below.
//            //        for (auto &elem: *f._range) {
//            //            _range->push_back(make_pair<>(copy(elem.first), copy(elem.second)));
//            //        }
//        }
//
//        for (auto &pair:*f._lterms) {
//            insert(pair.second);
//        }
//        for (auto &pair:*f._qterms) {
//            insert(pair.second);
//        }
//        for (auto &pair:*f._pterms) {
//            insert(pair.second);
//        }
//        if (f._val) {
//            _val = make_shared<vector<double>>(*f._val);
//        }
//        if (f._expr) {
//            if (f._expr->is_uexpr()) {
//                _expr = make_shared<uexpr>(*(uexpr*)(f._expr.get()));
//            }
//            else {
//                _expr = make_shared<bexpr>(*(bexpr*)(f._expr.get()));
//            }
//            embed(_expr);
//            //            _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//            _queue->push_back(_expr);
//        }
//        for (auto &df:*f._dfdx) {
//            (*_dfdx)[df.first] = make_shared<func_>(*df.second);
//        }
//        _nnz_j = f._nnz_j;
//        _nnz_h = f._nnz_h;
//        _hess_link = f._hess_link;
//        //        _nb_instances = f._nb_instances;
//        _nb_vars = f._nb_vars;
//        _indices = f._indices;
//        if (is_constant()) {
//            _evaluated = f._evaluated;
//        }
//
//        return *this;
//    }
//
//
//    func_& func_::operator=(func_&& f){
//        //        _evaluated = f._evaluated;
//        _to_str = f._to_str;
//
//        delete _all_range;
//        if (_vars) {
//            _vars->clear();
//            delete _vars;
//        }
//        if (_params){
//            _params->clear();
//            delete _params;
//        }
//
//        delete _range;
//        delete _convexity;
//        delete _sign;
//        delete _lterms;
//        delete _qterms;
//        delete _pterms;
//        delete _DAG;
//        delete _queue;
//        delete _cst;
//        set_type(func_c);
//        _ftype = f._ftype;
//        _return_type = f._return_type;
//        _all_convexity = f._all_convexity;
//        _all_sign = f._all_sign;
//        _all_range = f._all_range;
//        f._all_range = nullptr;
//        _lterms = f._lterms;
//        f._lterms = nullptr;
//        _qterms = f._qterms;
//        f._qterms = nullptr;
//        _pterms = f._pterms;
//        f._pterms = nullptr;
//        _expr = move(f._expr);
//        _DAG = f._DAG;
//        f._DAG = nullptr;
//        _queue = f._queue;
//        f._queue = nullptr;
//        _vars = f._vars;
//        f._vars = nullptr;
//        _params = f._params;
//        f._params = nullptr;
//        _cst = f._cst;
//        f._cst = nullptr;
//        _range = f._range;
//        f._range = nullptr;
//        _convexity = f._convexity;
//        f._convexity = nullptr;
//        _sign = f._sign;
//        f._sign = nullptr;
//        _is_transposed = f._is_transposed;
//        _is_vector = f._is_vector;
//        //        _is_matrix = f._is_matrix;
//        _is_constraint = f._is_constraint;
//        _is_hessian = f._is_hessian;
//        _dim[0] = f._dim[0];
//        _dim[1] = f._dim[1];
//        _embedded = f._embedded;
//        _nnz_j = f._nnz_j;
//        _nnz_h = f._nnz_h;
//        _hess_link = f._hess_link;
//        //        _nb_instances = f._nb_instances;
//        _nb_vars = f._nb_vars;
//        _dfdx = move(f._dfdx);
//        _val = move(f._val);
//        _indices = move(f._indices);
//        if (is_constant()) {
//            _evaluated = f._evaluated;
//        }
//        return *this;
//    }
//
//
//
//    func_::~func_(){
//        delete _all_range;
//        //        for (auto &pt: *_lterms) {
//        //            if (pt.second._coef->is_function()) {
//        //                delete pt.second._coef;
//        //            }
//        //        }
//        //        for (auto &pt: *_qterms) {
//        //            if (pt.second._coef->is_function()) {
//        //                delete pt.second._coef;
//        //            }
//        //        }
//        //        for (auto &pt: *_pterms) {
//        //            if (pt.second._coef->is_function()) {
//        //                delete pt.second._coef;
//        //            }
//        //        }
//        if (_vars) {
//            _vars->clear();
//            delete _vars;
//        }
//        if (_params) {
//            _params->clear();
//            delete _params;
//        }
//
//        delete _range;
//        delete _convexity;
//        delete _sign;
//        delete _lterms;
//        delete _qterms;
//        delete _pterms;
//        delete _DAG;
//        delete _queue;
//        delete _cst;
//    };
//
//    bool all_zeros(const string& s){
//        auto it = s.begin();
//        while (it != s.end()) {
//            if ((*it)!='-' && (*it)!='0' && (*it)!='.') {
//                return false;
//            }
//            it++;
//        }
//        return true;
//    }
//
//    bool constant_::is_zero() const{ /**< Returns true if constant equals 0 */
//        //        auto a = poly_to_str(this);
//        if (is_number() && all_zeros(poly_to_str(this))){
//            return true;
//        }
//        if (is_param() || is_var()) {
//            auto p_c = (param_*)this;
//            return p_c->get_all_sign()==zero_;
//        }
//        //        if (is_uexpr() || is_bexpr()) {
//        //            auto e_p = (expr*)this;
//        //            return e_p->get_all_sign()==zero_;
//        //        }
//        if (is_function()) {
//            return ((func_*)(this))->is_zero();
//        }
//        return false;
//    }
//
//    bool constant_::is_unit() const{ /**< Returns true if constant equals 1 */
//        if(is_number() && t_eval(this)==1 && !_is_vector && !_is_transposed){
//            return true;
//        }
//        //        if (is_param()) {
//        //            auto p_c = (param_*)this;
//        //            switch (p_c->get_intype()) {
//        //                case binary_:
//        //                    return ((param<bool>*)p_c)->is_unit();
//        //                    break;
//        //                case short_:
//        //                    return ((param<short>*)p_c)->is_unit();
//        //                    break;
//        //                case integer_:
//        //                    return ((param<int>*)p_c)->is_unit();
//        //                    break;
//        //                case float_:
//        //                    return ((param<float>*)p_c)->is_unit();
//        //                    break;
//        //                case double_:
//        //                    return ((param<double>*)p_c)->is_unit();
//        //                    break;
//        //                case long_:
//        //                    return ((param<long double>*)p_c)->is_unit();
//        //                    break;
//        //                default:
//        //                    break;
//        //            }
//        //        }
//        if (is_function()) {
//            return ((func_*)(this))->is_unit();
//        }
//        return false;
//    }
//
//
//    bool constant_::is_positive() const{
//        if (get_all_sign()==pos_) {
//            return true;
//        }
//        return false;
//    }
//
//
//    bool constant_::is_non_positive() const{
//        if (get_all_sign()==non_pos_) {
//            return true;
//        }
//        return false;
//    }
//
//    bool constant_::is_non_negative() const{
//        if (get_all_sign()==non_neg_) {
//            return true;
//        }
//        return false;
//    }
//
//    bool constant_::is_negative() const{
//        if (get_all_sign()==neg_) {
//            return true;
//        }
//        return false;
//    }
//
//    void func_::update_sign(){
//        //TODO
//    }
//
//
//    func_& func_::operator+=(const constant_& c){
//        if (c.is_zero()) {
//            return *this;
//        }
//        _evaluated = false;
//        if (c.is_number() || (!is_constant() && c.is_param())) {
//            _cst = add(_cst, c);
//            if (_cst->is_function()) {
//                embed(*(func_*)_cst);
//            }
//            update_sign(*_cst);
//            return *this;
//        }
//        if (c.is_param() || c.is_var()) {
//            this->insert(true, constant<double>(1), *(param_*)&c);
//        }
//        if (c.is_function()) {
//            func_* f = (func_*)&c;
//            if (is_constant() && !f->is_constant()) {//check _expr
//                return *this = (*f + *this);
//            }
//            if (!is_constant() && f->is_constant()) {//check _expr
//                _cst = add(_cst, c);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//                return *this;
//            }
//            if (!f->_cst->is_zero()) {
//                _cst = add(_cst, *f->_cst);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//            }
//            for (auto &pair:*f->_lterms) {
//                this->insert(pair.second);
//            }
//            for (auto &pair:*f->_qterms) {
//                this->insert(pair.second);
//            }
//            for (auto &pair:*f->_pterms) {
//                this->insert(pair.second);
//            }
//
//            if (_expr && f->_expr) {
//                _expr = make_shared<bexpr>(bexpr(plus_, make_shared<func_>((*_expr)), make_shared<func_>((*f->_expr))));
//                embed(_expr);
//                //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push_back(_expr);
//            }
//            else if (!_expr && f->_expr) {
//                if (f->_expr->is_uexpr()) {
//                    _expr = make_shared<uexpr>(*(uexpr*)(f->_expr.get()));
//                }
//                else {
//                    _expr = make_shared<bexpr>(*(bexpr*)(f->_expr.get()));
//                }
//                embed(_expr);
//                //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push_back(_expr);
//                if (!_vars->empty()) {
//                    _ftype = nlin_;
//                }
//            }
//            update_sign(*f);
//            if (_all_convexity==linear_ && f->_all_convexity==convex_) {
//                _all_convexity = convex_;
//            }
//            else if (_all_convexity==linear_ && f->_all_convexity==concave_) {
//                _all_convexity = concave_;
//            }
//            else if (_all_convexity==convex_ && f->_all_convexity==convex_) {
//                _all_convexity = convex_;
//            }
//            else if (_all_convexity==concave_ && f->_all_convexity==concave_) {
//                _all_convexity = concave_;
//            }
//            else if (_all_convexity==convex_ && (f->_all_convexity==convex_ || f->_all_convexity==undet_)) {
//                _all_convexity = undet_;
//            }
//            else if (_all_convexity==concave_ && (f->_all_convexity==concave_ || f->_all_convexity==undet_)) {
//                _all_convexity = undet_;
//            }
//            update_convexity();
//            return *this;
//        }
//        return *this;
//    }
//
//    func_& func_::operator-=(const constant_& c){
//        if (c.is_zero()) {
//            return *this;
//        }
//        _evaluated = false;
//        if (c.is_number() || (!is_constant() && c.is_param())) {
//            _cst = substract(_cst, c);
//            if (_cst->is_function()) {
//                embed(*(func_*)_cst);
//            }
//            update_sign(*_cst);
//            return *this;
//        }
//        if (c.is_param() || c.is_var()) {
//            this->insert(false, constant<double>(1), *(param_*)&c);
//        }
//        if (c.is_function()) {
//            func_* f = (func_*)&c;
//            if (!is_constant() && f->is_constant()) {
//                _cst = substract(_cst, c);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//                return *this;
//            }
//            if (is_constant() && !f->is_constant()) {//check _expr
//                return *this = (-1*(*f) + *this);
//            }
//            _cst = substract(_cst, *f->_cst);
//            if (_cst->is_function()) {
//                embed(*(func_*)_cst);
//            }
//            for (auto &pair:*f->_lterms) {
//                this->insert(!pair.second._sign, *pair.second._coef, *pair.second._p);
//            }
//            for (auto &pair:*f->_qterms) {
//                this->insert(!pair.second._sign, *pair.second._coef, *pair.second._p->first, *pair.second._p->second);
//            }
//            for (auto &pair:*f->_pterms) {
//                this->insert(!pair.second._sign, *pair.second._coef, *pair.second._l);
//            }
//            if (_expr && f->_expr) {
//                _expr = make_shared<bexpr>(bexpr(minus_, make_shared<func_>((*_expr)), make_shared<func_>((*f->_expr))));
//                embed(_expr);
//                //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push_back(_expr);
//            }
//            if (!_expr && f->_expr) {
//                if (f->_expr->is_uexpr()) {
//                    _expr = make_shared<uexpr>(*(uexpr*)(f->_expr.get()));
//                }
//                else {
//                    _expr = make_shared<bexpr>(*(bexpr*)(f->_expr.get()));
//                }
//                _expr->_coef *= -1;
//                embed(_expr);
//                //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push_back(_expr);
//                if (!_vars->empty()) {
//                    _ftype = nlin_;
//                }
//            }
//            update_sign(*f);
//            if (_all_convexity==linear_ && f->_all_convexity==concave_) {
//                _all_convexity = convex_;
//            }
//            else if (_all_convexity==linear_ && f->_all_convexity==convex_) {
//                _all_convexity = concave_;
//            }
//            else if (_all_convexity==convex_ && f->_all_convexity==concave_) {
//                _all_convexity = convex_;
//            }
//            else if (_all_convexity==concave_ && f->_all_convexity==convex_) {
//                _all_convexity = concave_;
//            }
//            else if (_all_convexity==convex_ && (f->_all_convexity==concave_ || f->_all_convexity==undet_)) {
//                _all_convexity = undet_;
//            }
//            else if (_all_convexity==concave_ && (f->_all_convexity==convex_ || f->_all_convexity==undet_)) {
//                _all_convexity = undet_;
//            }
//            update_convexity();
//            return *this;
//        }
//        return *this;
//    }
//
//
//    func_& func_::operator*=(const constant_& c){
//        if (is_zero()) {
//            return *this;
//        }
//        if (c.is_zero()) {
//            reset();
//            return *this;
//        }
//        if (is_unit()) {
//            *this = func_(c);
//            return *this;
//        }
//        if (c.is_unit()) {
//            return *this;
//        }
//
//        /* Case where c is a number */
//        if (c.is_number()){
//            if (!_cst->is_zero()) {
//                _cst = multiply(_cst, c);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//            }
//            for (auto &pair:*_lterms) {
//                pair.second._coef = multiply(pair.second._coef, c);
//            }
//            for (auto &pair:*_qterms) {
//                pair.second._coef = multiply(pair.second._coef, c);
//            }
//            for (auto &pair:*_pterms) {
//                pair.second._coef = multiply(pair.second._coef, c);
//            }
//            if (c.is_negative()) {
//                reverse_sign();
//            }
//            auto val = t_eval(&c);
//            if (_expr) {
//                _expr->_coef *= val;
//            }
//            if (_evaluated) {
//                for (unsigned i = 0; i<_val->size(); i++) {
//                    _val->at(i) *= val;
//                }
//            }
//            //            this->update_dot_dim(c);
//            return *this;
//        }
////        if (c.is_matrix() || is_matrix()) {
////            *this = func_(bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c)));
////            _evaluated = false;
////            return *this;
////        }
//        /* Case where the current function is not constant and the other operand is */
//        if(!is_constant() && (c.is_param() || (c.is_function() && ((func_*)&c)->is_constant()))) {
//            bool transp = false;
//            auto fc = func_(c);
//            if(is_linear() && _is_transposed){// Situation where var^T * c is transformed into (c^T*var)^T
//                fc.transpose();
//                this->transpose();
//                transp = true;
//            }
//            if (!_cst->is_zero()) {
//                _cst = multiply(_cst, fc);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//            }
//            for (auto &pair:*_lterms) {
//                pair.second._coef = multiply(pair.second._coef, fc);
//                if (pair.second._coef->is_function()) {
//                    embed(*(func_*)pair.second._coef);
//                }
//                //                update_nb_instances(pair.second);
//
//            }
//            for (auto &pair:*_qterms) {
//                pair.second._coef = multiply(pair.second._coef, fc);
//                if (pair.second._coef->is_function()) {
//                    embed(*(func_*)pair.second._coef);
//                }
//                //                update_nb_instances(pair.second);
//            }
//            for (auto &pair:*_pterms) {
//                pair.second._coef = multiply(pair.second._coef, fc);
//                if (pair.second._coef->is_function()) {
//                    embed(*(func_*)pair.second._coef);
//                }
//                //                update_nb_instances(pair.second); // TODO
//            }
//            if (_expr) {
//                if (c.is_function() && ((func_*)&c)->is_number()) {
//                    _expr->_coef *= t_eval(&c);
//                }
//                else {
//                    _expr = make_shared<bexpr>(bexpr(product_, make_shared<func_>((*_expr)), make_shared<func_>((c))));
//                    //                    _dim[0] = _expr->_dim[0];
//                    //                    _dim[1] = _expr->_dim[1];
//                    embed(_expr);
//                    //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                    _queue->push_back(_expr);
//                }
//            }
//            if (c.is_negative()) {
//                reverse_sign();
//            }
//            if (c.get_all_sign()==unknown_) {
//                _all_sign = unknown_;
//                if (!_qterms->empty()) {
//                    _all_convexity = undet_;
//                }
//            }
//            _evaluated = false;
//            if(transp){
//                this->transpose();
//            }
//            else {
//                this->update_dot(c);
//            }
//            return *this;
//        }
//        /* Case where the current function is constant and the other operand is not. */
//        if (is_constant() && (c.is_var() || (c.is_function() && !((func_*)&c)->is_constant()))) {
//            func_ f(c);
//            if (!f._cst->is_zero()) {
//                if (f._cst->is_function()) {
//                    auto fc = (func_*)f._cst;
//                    *fc = (*this)* (*fc);
//                    f.embed(*fc);
//                }
//                else {
//                    f._cst = multiply(f._cst, *this);
//                }
//            }
//            for (auto &pair:*f._lterms) {
//                if (pair.second._coef->is_function()) {
//                    auto fc = (func_*)pair.second._coef;
//                    *fc = (*this)* (*fc);
//                    f.embed(*fc);
//                }
//                else {
//                    pair.second._coef = multiply(pair.second._coef, *this);
//                }
//                if (pair.second._coef->_is_transposed) {
//                    pair.second._p->_is_vector = true;
//                    if (!pair.second._coef->is_number() && pair.second._coef->_dim[1]!=pair.second._p->_dim[0]) {
//                        DebugOn("vector dot product with mismatching dimensions, check your param/var dimensions");
//                    }
//                }
//                //                f.update_nb_instances(pair.second);
//
//            }
//            for (auto &pair:*f._qterms) {
//                if (pair.second._coef->is_function()) {
//                    auto fc = (func_*)pair.second._coef;
//                    *fc = (*this)* (*fc);
//                    f.embed(*fc);
//                }
//                else {
//                    pair.second._coef = multiply(pair.second._coef, *this);
//                }
//                if (pair.second._coef->_is_transposed) {
//                    pair.second._p->first->_is_vector = true;
//                    pair.second._p->second->_is_vector = true;
//                    if (!pair.second._coef->is_number() && pair.second._coef->_dim[1]!=pair.second._p->first->_dim[0]) {
//                        DebugOn("vector dot product with mismatching dimensions, check your param/var dimensions");
//                    }
//
//                }
//                //                update_nb_instances(pair.second);
//            }
//            for (auto &pair:*f._pterms) {
//                if (pair.second._coef->is_function()) {
//                    auto fc = (func_*)pair.second._coef;
//                    *fc = (*this)* (*fc);
//                    f.embed(*fc);
//                }
//                else {
//                    pair.second._coef = multiply(pair.second._coef, *this);
//                }
//                //                update_nb_instances(pair.second); // TODO
//            }
//            if (f._expr) {
//                if (this->is_number()) {
//                    f._expr->_coef *= t_eval(this);
//                }
//                else {
//                    f._expr = make_shared<bexpr>(bexpr(product_, make_shared<func_>(*this), make_shared<func_>((*f._expr))));
//                    f.embed(f._expr);
//                    //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                    f._queue->push_back(f._expr);
//                }
//            }
//            if (this->is_negative()) {
//                f.reverse_sign();
//            }
//            if (this->get_all_sign()==unknown_) {
//                f._all_sign = unknown_;
//                if (!f._qterms->empty()) {
//                    f._all_convexity = undet_;
//                }
//            }
//            f._evaluated = false;
//            if(update_dot(c)){
//                f._dim[0] = _dim[0];
//                f._dim[1] = _dim[1];
//            }
//            return *this = move(f);
//        }
//        if (c.is_param() || c.is_var()) {
////            if (c.is_matrix() || is_matrix()) {
////                *this = func_(bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c)));
////            }
////            else {
//            if(_is_transposed && !c._is_vector){
//                auto new_c = copy(c);
//                auto new_p = (param_*)new_c;
//                new_p->_is_vector = true;
//                new_p->_name = "["+new_p->_name+"]";
//                func_ f(*new_p);
//                *this *= f;
//                delete new_c;
//            }
//            else {
//                func_ f(c);
//                *this *= f;
//            }
////            }
//            _evaluated = false;
//            return *this;
//        }
//        if (_expr || (c.is_function() && ((func_*)&c)->_expr)) {
//            auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//            *this = func_(be);
//            _evaluated = false;
//            return *this;
//        }
//        /* Case where the multiplication invlolves multiplying variables/parameters together, i.e., they are both parametric or both include variables  */
//        if (c.is_function()) {
//            func_* f = (func_*)&c;
//            constant_* coef;
//            vector<bool>* is_sum = nullptr;
//            func_ res;
//            for (auto& t1: *_pterms) {
//                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), we cannot factor the coefficients. Just create a binary expression and return it.
//                    auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                    *this = func_(be);
//                    _evaluated = false;
//                    return *this;
//                }
//                for (auto& t2: *f->_pterms) {
//                    is_sum = nullptr;
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    auto newl(*t1.second._l);
//                    for (auto& it: *t2.second._l) {
//                        newl.push_back(make_pair<>(it.first, it.second));
//                    }
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                for (auto& t2: *f->_qterms) {
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function), see comment above.
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    auto newl(*t1.second._l);
//                    newl.push_back(make_pair<>((t2.second._p->first), 1));
//                    newl.push_back(make_pair<>((t2.second._p->second), 1));
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                for (auto& t2: *f->_lterms) {
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial function)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    auto newl(*t1.second._l);
//                    newl.push_back(make_pair<>((t2.second._p), 1));
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                if (!f->_cst->is_zero()) {
//                    auto newl(*t1.second._l);
//                    coef = copy(*f->_cst);
//                    coef = multiply(coef, *t1.second._coef);
//                    res.insert(t1.second._sign, *coef, newl);
//                    delete coef;
//                }
//            }
//
//            for (auto& t1: *_qterms) {
//                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
//                    auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                    *this = func_(be);
//                    _evaluated = false;
//                    return *this;
//                }
//                for (auto& t2: *f->_pterms) {
//                    is_sum = nullptr;
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    auto newl(*t2.second._l);
//                    newl.push_front(make_pair<>(t1.second._p->first, 1));
//                    newl.push_front(make_pair<>(t1.second._p->second, 1));
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                for (auto& t2: *f->_qterms) {
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    list<pair<param_*, int>> newl;
//                    newl.push_back(make_pair<>(t1.second._p->first, 1));
//                    newl.push_back(make_pair<>(t1.second._p->second, 1));
//                    newl.push_back(make_pair<>(t2.second._p->first, 1));
//                    newl.push_back(make_pair<>(t2.second._p->second, 1));
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                for (auto& t2: *f->_lterms) {
//                    is_sum = nullptr;
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    list<pair<param_*, int>> newl;
//                    newl.push_back(make_pair<>(t1.second._p->first, 1));
//                    newl.push_back(make_pair<>(t1.second._p->second, 1));
//                    newl.push_back(make_pair<>(t2.second._p, 1));
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                if (!f->_cst->is_zero()) {
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *f->_cst);
//                    res.insert(t1.second._sign, *coef, *t1.second._p->first, *t1.second._p->second);
//                    delete coef;
//                }
//
//            }
//            for (auto& t1: *_lterms) {
////                if (t1.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(Quadratic term)
////                    auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
////                    *this = func_(be);
////                    _evaluated = false;
////                    return *this;
////                }
//                for (auto& t2: *f->_pterms) {
//                    is_sum = nullptr;
//                    if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    auto newl(*t2.second._l);
//                    newl.push_front(make_pair<>((t1.second._p), 1));
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                for (auto& t2: *f->_qterms) {
//                    if (t1.second._coef->_is_transposed || t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    list<pair<param_*, int>> newl;
//                    newl.push_back(make_pair<>(t1.second._p, 1));
//                    newl.push_back(make_pair<>(t2.second._p->first, 1));
//                    newl.push_back(make_pair<>(t2.second._p->second, 1));
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, newl);
//                    delete coef;
//                }
//                for (auto& t2: *f->_lterms) {
////                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
////                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
////                        *this = func_(be);
////                        _evaluated = false;
////                        return *this;
////                    }
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(!(t1.second._sign^t2.second._sign), *coef, *t1.second._p, *t2.second._p, _is_transposed);
//                    delete coef;
//                }
//                if (!f->_cst->is_zero()) {
//                    coef = copy(*t1.second._coef);
//                    coef = multiply(coef, *f->_cst);
//                    res.insert(t1.second._sign, *coef, *t1.second._p);
//                    delete coef;
//                }
//            }
//            if (!_cst->is_zero()) {
//                for (auto& t2: *f->_pterms) {
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    coef = copy(*_cst);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(t2.second._sign, *coef, *t2.second._l);
//                    delete coef;
//                }
//                for (auto& t2: *f->_qterms) {
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    coef = copy(*_cst);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(t2.second._sign, *coef, *t2.second._p->first, *t2.second._p->second);
//                    delete coef;
//                }
//                for (auto& t2: *f->_lterms) {
//                    if (t2.second._coef->_is_transposed) {// If the coefficient in front is transposed: a^T.(polynomial term)
//                        auto be = bexpr(product_, make_shared<func_>(*this), make_shared<func_>(c));
//                        *this = func_(be);
//                        _evaluated = false;
//                        return *this;
//                    }
//                    coef = copy(*_cst);
//                    coef = multiply(coef, *t2.second._coef);
//                    res.insert(t2.second._sign, *coef, *t2.second._p);
//                    delete coef;
//                }
//                if (!f->_cst->is_zero()) {
//                    coef = copy(*_cst);
//                    coef = multiply(coef, *f->_cst);
//                    delete res._cst;
//                    res._cst = coef;
//                    if (_cst->is_function()) {
//                        embed(*(func_*)_cst);
//                    }
//                }
//            }
//            res.update_dot_dim(*this, c);
//            *this = move(res);
//        }
//        _evaluated = false;
//        return *this;
//    }
//
//    func_& func_::operator/=(const constant_& c){
//        if (is_zero()) {
//            return *this;
//        }
//        if (c.is_unit()) {
//            return *this;
//        }
//        if (c.is_zero()) {
//            throw invalid_argument("dividing by zero!\n");
//        }
//        _evaluated = false;
//        /* Case where c is a number */
//        if (c.is_number()){
//            if (!_cst->is_zero()) {
//                _cst = divide(_cst, c);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//            }
//            for (auto &pair:*_lterms) {
//                pair.second._coef = divide(pair.second._coef, c);
//            }
//            for (auto &pair:*_qterms) {
//                pair.second._coef = divide(pair.second._coef, c);
//            }
//            for (auto &pair:*_pterms) {
//                pair.second._coef = divide(pair.second._coef, c);
//            }
//            if (_expr) {
//                _expr = make_shared<bexpr>(bexpr(div_, make_shared<func_>((*_expr)), make_shared<func_>((c))));
//                embed(_expr);
//                //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push_back(_expr);
//            }
//            if (c.is_negative()) {
//                reverse_convexity();
//                reverse_sign();
//            }
//        }
//        /* Case where the current function is not constant and the other operand is */
//        if(!is_constant() && (c.is_param() || (c.is_function() && ((func_*)&c)->is_constant()))) {
//            if (!_cst->is_zero()) {
//                _cst = divide(_cst, c);
//                if (_cst->is_function()) {
//                    embed(*(func_*)_cst);
//                }
//            }
//            for (auto &pair:*_lterms) {
//                pair.second._coef = divide(pair.second._coef, c);
//            }
//            for (auto &pair:*_qterms) {
//                pair.second._coef = divide(pair.second._coef, c);
//            }
//            for (auto &pair:*_pterms) {
//                pair.second._coef = divide(pair.second._coef, c);
//            }
//            if (c.is_negative()) {
//                reverse_convexity();
//                reverse_sign();
//            }
//            if (c.get_all_sign()==unknown_) {
//                _all_sign = unknown_;
//                if (!_qterms->empty()) {
//                    _all_convexity = undet_;
//                }
//            }
//            if (_expr) {
//                _expr = make_shared<bexpr>(bexpr(div_, make_shared<func_>((*_expr)), make_shared<func_>((c))));
//                embed(_expr);
//                //                _DAG->insert(make_pair<>(_expr->get_str(), _expr));
//                _queue->push_back(_expr);
//            }
//            return *this;
//        }
//        auto be = bexpr(div_, make_shared<func_>(*this), make_shared<func_>((c)));
//        *this = func_(be);
//        return *this;
//    }
//
//    //    void func_::update_nb_ind(){
//    //        _nb_instances = 0;
//    //        for (auto &p: *_vars) {
//    //            if (!p.second.first->_is_vector) {
//    //                _nb_instances = max(_nb_instances, p.second.first->get_dim());
//    //            }
//    //        }
//    //    }
//
//    
//
    

    
    
    
    
    /**
     Reverse the convexity property of the current function
     */
    void func_::reverse_convexity(){
        if (_all_convexity==convex_) {
            _all_convexity=concave_;
        }
        else if (_all_convexity==concave_) {
            _all_convexity=convex_;
        }
    }

    

//
    void func_::update_sign_add(const constant_& c){
        _all_sign = sign_add(_all_sign, c.get_all_sign());
    }
    
    void func_::update_sign_multiply(const constant_& c){
        _all_sign = sign_product(_all_sign, c.get_all_sign());
    }
//
//    void func_::set_max_dim(const lterm& l){
//        //        assert(_dim[0] <= l._p->_dim[0] && _dim[1] <= l._p->_dim[1]);
//        if(l._p->is_indexed() && l._p->_indices->_ids->size()>1){//This is a multi-instance product
//            return;
//        }
//        if(!l._coef->_is_transposed && !l._p->_is_vector){
//            _dim[0] = max(_dim[0], l._p->_dim[0]);
//        }
//        if(l._coef->_is_vector){
//            _dim[0] = l._coef->_dim[0];
//        }
//        if(l._p->_is_vector){
//            _dim[1] = l._p->_dim[1];
//        }
//        /* Instructions above work for matrix and vector dot products, below we check if it's a component-wise vector,matrix product */
//        if(!l._coef->is_matrix() && l._p->is_matrix()){
//            _dim[0] = l._p->_dim[0];
//        }
//        if(l._coef->is_matrix() && !l._p->is_matrix() && l._p->_is_transposed){
//            _dim[1] = l._coef->_dim[1];
//        }
//    }
//
//    void func_::set_max_dim(const qterm& q){
//        if(q._p->first->is_indexed() && q._p->first->_indices->_ids->size()>1){//This is a multi-instance product
//            return;
//        }
//        if(q._c_p1_transposed){ // q = (coef*p1)^T*p2
//            _dim[0] = max(_dim[0],q._p->first->_dim[1]);
//            _dim[1] = max(_dim[1],q._p->second->_dim[1]);
//            return;
//        }
//        if(!q._coef->_is_transposed && !q._p->first->_is_vector){
//            _dim[0] = max(_dim[0], q._p->first->_dim[0]);
//        }
//        if (q._p->first->_is_vector) {
//            _dim[0] = q._p->first->_dim[0];
//        }
//        if (!q._coef->_is_transposed && q._p->second->_is_vector) {
//            _dim[1] = q._p->second->_dim[1];
//        }
//        /* Instructions above work for matrix and vector dot products, below we check if it's a component-wise vector,matrix product */
//        if(!q._p->first->is_matrix() && q._p->second->is_matrix()){
//            _dim[0] = q._p->second->_dim[0];
//        }
//        if(q._p->first->is_matrix() && !q._p->second->is_matrix() && q._p->second->_is_transposed){
//            _dim[1] = q._p->first->_dim[1];
//        }
//    }
//
//
//    void func_::update_sign(const lterm& l){
//        Sign sign = get_all_sign(l);
//        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
//            _all_sign = unknown_;
//        }
//        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
//            _all_sign = unknown_;
//        }
//        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
//            _all_sign = sign;
//        }
//    }
//
//    void func_::update_sign(const qterm& q){
//        Sign sign = get_all_sign(q);
//        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
//            _all_sign = unknown_;
//        }
//        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
//            _all_sign = unknown_;
//        }
//        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
//            _all_sign = sign;
//        }
//    }
//
//    void func_::update_sign(const pterm& p){
//        Sign sign = get_all_sign(p);
//        if (sign==unknown_ || ((_all_sign==non_neg_ || _all_sign==pos_) && sign!=non_neg_ && sign!=pos_)) {
//            _all_sign = unknown_;
//        }
//        else if((_all_sign==non_pos_ || _all_sign==neg_) && sign!=non_pos_ && sign!=neg_){
//            _all_sign = unknown_;
//        }
//        else if(_all_sign==zero_ || _all_sign==pos_ || _all_sign==neg_){// take weaker sign
//            _all_sign = sign;
//        }
//    }
//
//    void func_::update_convexity(const qterm& q){
//        Convexity conv = get_convexity(q);
//        if (_all_convexity==undet_ || conv ==undet_ || (_all_convexity==convex_ && conv==concave_) || (_all_convexity==concave_ && conv==convex_)) {
//            _all_convexity = undet_;
//        }
//        else {
//            _all_convexity = conv;
//        }
//    }
//
//    bool func_::insert(bool sign, const constant_& coef, const param_& p);

    
//
//    bool func_::insert(bool sign, const constant_& coef, const param_& p1, const param_& p2, bool c_p1_transposed);

//    bool func_::insert(const qterm& term)

//    bool func_::insert(bool sign, const constant_& coef, const param_& p, int exp)

//    bool func_::insert(bool sign, const constant_& coef, const list<pair<shared_ptr<param_>, int>>& l)

    
//
//    //    void func_::insert(expr& e){
//    //    //    insert(term._sign, *term._coef, *term._l);
//    //        auto name = e.get_str();
//    ////        auto pair_it = _DAG->find(name);
//    //
//    //        _ftype = nlin_;
//    //
//    //    //    if (pair_it == _DAG->end()) {
//    //    //        update_sign(e);
//    //    //        update_convexity(e);
//    //    //        _DAG->insert(make_pair<>(name, e));
//    //    //    }
//    //    //    else {
//    //    //        if (pair_it->second._sign == e.sign) {
//    //    ////            pair_it->second._coef = add(pair_it->second._coef, coef);
//    //    ////            if (!pair_it->second._sign) { // both negative
//    //    ////                pair_it->second._sign = true;
//    //    ////            }
//    //    //        }
//    //    //        else{
//    //    ////            pair_it->second._coef = substract(pair_it->second._coef, coef);
//    //    //        }
//    //    ////        if (pair_it->second._coef->is_zero()) {
//    //    ////            if (p1.is_var()) {
//    //    ////                decr_occ_var(s1);
//    //    ////            }
//    //    ////            else {
//    //    ////                decr_occ_param(s1);
//    //    ////            }
//    //    ////            if (p2.is_var()) {
//    //    ////                decr_occ_var(s2);
//    //    ////            }
//    //    ////            else {
//    //    ////                decr_occ_param(s2);
//    //    ////            }
//    //    ////            _qterms->erase(pair_it);
//    //    ////            update_sign();
//    //    ////            update_convexity();
//    //    ////        }
//    //    ////        else {
//    //    ////            update_sign(pair_it->second);
//    //    ////            update_convexity(pair_it->second);
//    //    ////        }
//    //    ////        return false;
//    //    //    }
//    //    ////    _DAG->insert(<#const value_type &__v#>)
//    //    }
//
//
//
//
//    func_ cos(const constant_& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::cos(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        res._expr = make_shared<uexpr>(uexpr(cos_, make_shared<func_>((c))));
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//    func_ cos(constant_&& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::cos(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        res._expr = make_shared<uexpr>(uexpr(cos_, make_shared<func_>((move(c)))));
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//
//    func_ sin(const constant_& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::sin(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        res._expr = make_shared<uexpr>(uexpr(sin_, make_shared<func_>((c))));
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//    func_ sin(constant_&& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::sin(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        res._expr = make_shared<uexpr>(uexpr(sin_, make_shared<func_>((move(c)))));
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//
//    func_ sqrt(const constant_& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::sqrt(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        auto exp = make_shared<uexpr>(uexpr(sqrt_,make_shared<func_>((c))));
//        res._expr = exp;
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//            if (exp->_son->is_concave()) {
//                res._all_convexity = concave_;
//            }
//            else {
//                res._all_convexity = undet_;
//            }
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//    func_ sqrt(constant_&& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::sqrt(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        auto exp = make_shared<uexpr>(uexpr(sqrt_, make_shared<func_>((move(c)))));
//        res._expr = exp;
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//            if (exp->_son->is_concave()) {
//                res._all_convexity = concave_;
//            }
//            else {
//                res._all_convexity = undet_;
//            }
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//
//    func_ expo(const constant_& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::exp(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        res._is_vector = c._is_vector;
//        //        res._is_matrix = c._is_matrix;
//        res._is_transposed = c._is_transposed;
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        auto exp = make_shared<uexpr>(uexpr(exp_, make_shared<func_>((c))));
//        res._expr = exp;
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//            if (exp->_son->is_convex()) {
//                res._all_convexity = convex_;
//            }
//            else {
//                res._all_convexity = undet_;
//            }
//        }
//        return res;
//    };
//
//
//    func_ expo(constant_&& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::exp(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        res._is_vector = c._is_vector;
//        //        res._is_matrix = c._is_matrix;
//        res._is_transposed = c._is_transposed;
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        auto exp = make_shared<uexpr>(uexpr(exp_, make_shared<func_>((move(c)))));
//        res._expr = exp;
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//            if (exp->_son->is_convex()) {
//                res._all_convexity = convex_;
//            }
//            else {
//                res._all_convexity = undet_;
//            }
//        }
//        return res;
//    };
//
//
//
//    func_ log(const constant_& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::log(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//
//
//        auto exp = make_shared<uexpr>(uexpr(log_, make_shared<func_>((c))));
//        res._expr = exp;
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//            if (exp->_son->is_concave()) {
//                res._all_convexity = concave_;
//            }
//            else {
//                res._all_convexity = undet_;
//            }
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//    func_ log(constant_&& c){
//        func_ res;
//        if (c.is_number()) {
//            res._val->resize(1, std::log(t_eval(&c)));
//        }
//        else {
//            res._evaluated = false;
//        }
//        auto exp = make_shared<uexpr>(uexpr(log_, make_shared<func_>((move(c)))));
//        res._expr = exp;
//        res.embed(res._expr);
//        //        res._DAG->insert(make_pair<>(res._expr->get_str(), res._expr));
//        res._queue->push_back(res._expr);
//        if (!res._vars->empty()) {
//            res._ftype = nlin_;
//            if (exp->_son->is_concave()) {
//                res._all_convexity = concave_;
//            }
//            else {
//                res._all_convexity = undet_;
//            }
//        }
//        res._dim[0] = c._dim[0];
//        res._dim[1] = c._dim[1];
//        return res;
//    };
//
//    constant_* copy(const constant_& c2){/**< Return a copy c2 detecting the right class, i.e., constant<>, param<>, uexpr , bexpr or func. WARNING: allocates memory! */
//
//        switch (c2.get_type()) {
//            case binary_c: {
//                return new constant<bool>(*(constant<bool>*)(&c2));
//                break;
//            }
//            case short_c: {
//                return new constant<short>(*(constant<short>*)(&c2));
//                break;
//            }
//            case integer_c: {
//                return new constant<int>(*(constant<int>*)(&c2));
//                break;
//            }
//            case float_c: {
//                return new constant<float>(*(constant<float>*)(&c2));
//                break;
//            }
//            case double_c: {
//                return new constant<double>(*(constant<double>*)(&c2));
//                break;
//            }
//            case long_c: {
//                return new constant<long double>(*(constant<long double>*)(&c2));
//                break;
//            }
//            case complex_c: {
//                return new constant<Cpx>(*(constant<Cpx>*)(&c2));
//                break;
//            }
//            case par_c:{
//                auto p_c2 = (param_*)(&c2);
//                switch (p_c2->get_intype()) {
//                    case binary_:
//                        return new param<bool>(*(param<bool>*)p_c2);
//                        break;
//                    case short_:
//                        return new param<short>(*(param<short>*)p_c2);
//                        break;
//                    case integer_:
//                        return new param<int>(*(param<int>*)p_c2);
//                        break;
//                    case float_:
//                        return new param<float>(*(param<float>*)p_c2);
//                        break;
//                    case double_:
//                        return new param<double>(*(param<double>*)p_c2);
//                        break;
//                    case long_:
//                        return new param<long double>(*(param<long double>*)p_c2);
//                        break;
//                    case complex_:
//                        return new param<Cpx>(*(param<Cpx>*)p_c2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case uexp_c: {
//                return new uexpr(*(uexpr*)&c2);
//                break;
//            }
//            case bexp_c: {
//                return new bexpr(*(bexpr*)&c2);
//                break;
//            }
//            case var_c:{
//                auto p_c2 = (param_*)(&c2);
//                switch (p_c2->get_intype()) {
//                    case binary_:
//                        return new var<bool>(*(var<bool>*)p_c2);
//                        break;
//                    case short_:
//                        return new var<short>(*(var<short>*)p_c2);
//                        break;
//                    case integer_:
//                        return new var<int>(*(var<int>*)p_c2);
//                        break;
//                    case float_:
//                        return new var<float>(*(var<float>*)p_c2);
//                        break;
//                    case double_:
//                        return new var<double>(*(var<double>*)p_c2);
//                        break;
//                    case long_:
//                        return new var<long double>(*(var<long double>*)p_c2);
//                        break;
//                    case complex_:
//                        return new var<Cpx>(*(var<Cpx>*)p_c2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case func_c: {
//                return new func_(c2);
//                break;
//            }
//
//            default:
//                break;
//        }
//        return nullptr;
//    }
//
//    bool equals(const constant_* c1, const constant_* c2){/**< Checks if c2 equals c1 detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//        if ((!c1 && !c2) || (c1==c2)) {
//            return true;
//        }
//        if ((c1 && !c2) || (!c1 && c2)) {
//            return false;
//        }
//
//        switch (c2->get_type()) {
//            case binary_c: {
//                return (c1->is_binary() && *(constant<bool>*)c1 == *(constant<bool>*)c2);
//                break;
//            }
//            case short_c: {
//                return (c1->is_short() && *(constant<short>*)c1 == *(constant<short>*)c2);
//                break;
//            }
//            case integer_c: {
//                return (c1->is_integer() && *(constant<int>*)c1 == *(constant<int>*)c2);
//                break;
//            }
//            case float_c: {
//                return (c1->is_float() && *(constant<float>*)c1 == *(constant<float>*)c2);
//                break;
//            }
//            case double_c: {
//                return (c1->is_double() && *(constant<double>*)c1 == *(constant<double>*)c2);
//                break;
//            }
//            case long_c: {
//                return (c1->is_long() && *(constant<long double>*)c1 == *(constant<long double>*)c2);
//                break;
//            }
//            case complex_c: {
//                return (c1->is_complex() && *(constant<Cpx>*)c1 == *(constant<Cpx>*)c2);
//                break;
//            }
//            case par_c:{
//                return (c1->is_param() && *(param_ *)c1 == *(param_ *)c2);
//                break;
//            }
//            case uexp_c: {
//                return (c1->is_uexpr() && *(uexpr *)c1 == *(uexpr *)c2);
//                break;
//            }
//            case bexp_c: {
//                return (c1->is_bexpr() && *(bexpr *)c1 == *(bexpr *)c2);
//                break;
//            }
//            case var_c:{
//                return (c1->is_var() && *(param_ *)c1 == *(param_ *)c2);
//                break;
//            }
//            case func_c:{
//                if (c1->is_function()){
//                    auto f1 = (func_ *)c1;
//                    auto f2 = (func_ *)c2;
//                    return *f1==*f2;
//                }
//                break;
//            }
//            default:
//                break;
//        }
//        return false;
//    }
//
//    void set_val(param_* p, double val, size_t i){ /**< Polymorphic set_val */
//        switch (p->get_intype()) {
//            case binary_:
//                ((param<bool>*)p)->set_val(i, val);
//                break;
//            case short_:
//                ((param<short>*)p)->set_val(i, val);
//                break;
//            case integer_:
//                ((param<int>*)p)->set_val(i, val);
//                break;
//            case float_:
//                ((param<float>*)p)->set_val(i, val);
//                break;
//            case double_:
//                ((param<double>*)p)->set_val(i, val);
//                break;
//            case long_:
//                ((param<long double>*)p)->set_val(i, val);
//                break;
//            default:
//                throw invalid_argument("Polymorphic set_val only works for arithmetic types");
//                break;
//        }
//    }
//
//
//    /* Polymorphic functions */
//
//    //void reverse_sign(constant_* c){ /**< Reverses the sign of the constant. */
//    //    switch (c->get_type()) {
//    //        case binary_c: {
//    //            ((constant<bool>*)c)->set_val(!((constant<bool>*)c)->eval());
//    //            break;
//    //        }
//    //        case short_c: {
//    //            ((constant<short>*)c)->set_val(-1*((constant<short>*)c)->eval());
//    //            break;
//    //        }
//    //        case integer_c: {
//    //            ((constant<int>*)c)->set_val(-1*((constant<int>*)c)->eval());
//    //            break;
//    //        }
//    //        case float_c: {
//    //            ((constant<float>*)c)->set_val(-1*((constant<float>*)c)->eval());
//    //            break;
//    //        }
//    //        case double_c: {
//    //            ((constant<double>*)c)->set_val(-1*((constant<double>*)c)->eval());
//    //            break;
//    //        }
//    //        case long_c: {
//    //            ((constant<long double>*)c)->set_val(-1*((constant<long double>*)c)->eval());
//    //            break;
//    //        }
//    //        default:
//    //            throw invalid_argument("Cannot reverse sign of non-numeric constant");
//    //            break;
//    //    }
//    //
//    //}
//
//    void poly_print(const constant_* c){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//        cout << poly_to_str(c);
//    }
//
//    string matrix_to_str(const constant_* c){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//
//        if (!c) {
//            return "null";
//        }
//        if(!c->is_matrix()){
//            return poly_to_str(c);
//        }
//        switch (c->get_type()) {
//            case par_c:{
//                return ((param_*)(c))->get_name(false,false);
//                break;
//            }
//            case uexp_c: {
//                return ((uexpr*)c)->get_str();
//                break;
//            }
//            case bexp_c: {
//                return ((bexpr*)c)->get_str();
//                break;
//            }
//            case var_c: {
//                return ((param_*)(c))->get_name(false,false);
//                break;
//            }
//            case func_c: {
//                return ((func_*)c)->to_str();
//                break;
//            }
//            default:
//                break;
//        }
//        return "null";
//    }
//
//    string poly_to_str(const constant_* c){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//
//        if (!c) {
//            return "null";
//        }
//        switch (c->get_type()) {
//            case binary_c: {
//                return ((constant<bool>*)(c))->to_str();
//                break;
//            }
//            case short_c: {
//                return ((constant<short>*)(c))->to_str();
//                break;
//            }
//            case integer_c: {
//                return ((constant<int>*)(c))->to_str();
//                break;
//            }
//            case float_c: {
//                return ((constant<float>*)(c))->to_str();
//                break;
//            }
//            case double_c: {
//                return ((constant<double>*)(c))->to_str();
//                break;
//            }
//            case long_c: {
//                return ((constant<long double>*)(c))->to_str();
//                break;
//            }
//            case complex_c: {
//                return ((constant<Cpx>*)(c))->to_str();
//                break;
//            }
//            case par_c:{
//                return ((param_*)(c))->get_name(false,false);
//                break;
//            }
//            case uexp_c: {
//                return ((uexpr*)c)->get_str();
//                break;
//            }
//            case bexp_c: {
//                return ((bexpr*)c)->get_str();
//                break;
//            }
//            case var_c: {
//                return ((param_*)(c))->get_name(false,false);
//                break;
//            }
//            case func_c: {
//                return ((func_*)c)->to_str();
//                break;
//            }
//            default:
//                break;
//        }
//        return "null";
//    }
//
//    //    string poly_to_str(const constant_* c, size_t inst1, size_t inst2){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//    //
//    //        if (!c) {
//    //            return "null";
//    //        }
//    //        switch (c->get_type()) {
//    //            case binary_c: {
//    //                return ((constant<bool>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case short_c: {
//    //                return ((constant<short>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case integer_c: {
//    //                return ((constant<int>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case float_c: {
//    //                return ((constant<float>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case double_c: {
//    //                return ((constant<double>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case long_c: {
//    //                return ((constant<long double>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case complex_c: {
//    //                return ((constant<Cpx>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case par_c:{
//    //                auto p_c = (param_*)(c);
//    //                switch (p_c->get_intype()) {
//    //                    case binary_:
//    //                        return ((param<bool>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    case short_:
//    //                        return ((param<short>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    case integer_:
//    //                        return ((param<int>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    case float_:
//    //                        return ((param<float>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    case double_:
//    //                        return ((param<double>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    case long_:
//    //                        return ((param<long double>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    case complex_:
//    //                        return ((param<Cpx>*)p_c)->to_str(inst1,inst2);
//    //                        break;
//    //                    default:
//    //                        break;
//    //                }
//    //                break;
//    //            }
//    //            case uexp_c: {
//    //                return ((uexpr*)c)->to_str(inst2);
//    //                break;
//    //            }
//    //            case bexp_c: {
//    //                return ((bexpr*)c)->to_str(inst2);
//    //                break;
//    //            }
//    //            case var_c: {
//    //                auto p_c = (param_*)(c);
//    //                switch (p_c->get_intype()) {
//    //                    case binary_:
//    //                        return ((var<bool>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    case short_:
//    //                        return ((var<short>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    case integer_:
//    //                        return ((var<int>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    case float_:
//    //                        return ((var<float>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    case double_:
//    //                        return ((var<double>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    case long_:
//    //                        return ((var<long double>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    case complex_:
//    //                        return ((var<Cpx>*)p_c)->get_name(inst1,inst2);
//    //                        break;
//    //                    default:
//    //                        break;
//    //                }
//    //                break;
//    //            }
//    //
//    //            case func_c: {
//    //                return ((func_*)c)->to_str(inst2);
//    //                break;
//    //            }
//    //            default:
//    //                break;
//    //        }
//    //        return "null";
//    //    }
//    //
//    //
//    //    string poly_to_str(const constant_* c, size_t inst){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//    //
//    //        if (!c) {
//    //            return "null";
//    //        }
//    //        switch (c->get_type()) {
//    //            case binary_c: {
//    //                return ((constant<bool>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case short_c: {
//    //                return ((constant<short>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case integer_c: {
//    //                return ((constant<int>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case float_c: {
//    //                return ((constant<float>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case double_c: {
//    //                return ((constant<double>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case long_c: {
//    //                return ((constant<long double>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case complex_c: {
//    //                return ((constant<Cpx>*)(c))->to_str();
//    //                break;
//    //            }
//    //            case par_c:{
//    //                auto p_c = (param_*)(c);
//    //                switch (p_c->get_intype()) {
//    //                    case binary_:
//    //                        return ((param<bool>*)p_c)->to_str(inst);
//    //                        break;
//    //                    case short_:
//    //                        return ((param<short>*)p_c)->to_str(inst);
//    //                        break;
//    //                    case integer_:
//    //                        return ((param<int>*)p_c)->to_str(inst);
//    //                        break;
//    //                    case float_:
//    //                        return ((param<float>*)p_c)->to_str(inst);
//    //                        break;
//    //                    case double_:
//    //                        return ((param<double>*)p_c)->to_str(inst);
//    //                        break;
//    //                    case long_:
//    //                        return ((param<long double>*)p_c)->to_str(inst);
//    //                        break;
//    //                    case complex_:
//    //                        return ((param<Cpx>*)p_c)->to_str(inst);
//    //                        break;
//    //                    default:
//    //                        break;
//    //                }
//    //                break;
//    //            }
//    //            case uexp_c: {
//    //                return ((uexpr*)c)->to_str(inst);
//    //                break;
//    //            }
//    //            case bexp_c: {
//    //                return ((bexpr*)c)->to_str(inst);
//    //                break;
//    //            }
//    //            case var_c: {
//    //                auto p_c = (param_*)(c);
//    //                switch (p_c->get_intype()) {
//    //                    case binary_:
//    //                        return ((var<bool>*)p_c)->get_name(inst);
//    //                        break;
//    //                    case short_:
//    //                        return ((var<short>*)p_c)->get_name(inst);
//    //                        break;
//    //                    case integer_:
//    //                        return ((var<int>*)p_c)->get_name(inst);
//    //                        break;
//    //                    case float_:
//    //                        return ((var<float>*)p_c)->get_name(inst);
//    //                        break;
//    //                    case double_:
//    //                        return ((var<double>*)p_c)->get_name(inst);
//    //                        break;
//    //                    case long_:
//    //                        return ((var<long double>*)p_c)->get_name(inst);
//    //                        break;
//    //                    case complex_:
//    //                        return ((var<Cpx>*)p_c)->get_name(inst);
//    //                        break;
//    //                    default:
//    //                        break;
//    //                }
//    //                break;
//    //            }
//    //            case func_c: {
//    //                return ((func_*)c)->to_str(inst);
//    //                break;
//    //            }
//    //            default:
//    //                break;
//    //        }
//    //        return "null";
//    //    }
//
//
//
//
//    func_ operator+(const constant_& c1, const constant_& c2){
//        return func_(c1) += c2;
//    }
//
//    //    func_ operator+(func_&& f, const constant_& c){
//    //        return f += c;
//    //    }
//    //
//    //    func_ operator+(const constant_& c, func_&& f){
//    //        return f += c;
//    //    }
//
//
//    func_ operator-(const constant_& c1, const constant_& c2){
//        return func_(c1) -= c2;
//    }
//
//    //    func_ operator-(func_&& f, const constant_& c){
//    //        return f -= c;
//    //    }
//    //
//    //    func_ operator-(const constant_& c, func_&& f){
//    //        return (f *= -1) += c;
//    //    }
//    //
//
//    func_ operator*(const constant_& c1, const constant_& c2){// Rewrite this to change res after the multiplication is done, make sure both vars are now vecs.
//        if(c1.is_number()) {
//            if (c1._is_transposed && !c1.is_matrix()) {//TODO check when matrix
//                auto new_c2 = copy(c2);
//                new_c2->_is_vector = true;
//                auto res = func_(c1);
//                res._dim[1] = c2._dim[0];
//                res *= move(*new_c2);
//                res._is_vector = false;
//                res._is_transposed = false;
//                delete new_c2;
//                return res;
//            }
//            else {
//                return func_(c1) *= c2;
//            }
//        }
//        else {
////        if (c1._is_transposed) {//TODO check when matrix
////            auto new_c2 = copy(c2);
////            new_c2->_is_vector = true;
////            return func_(c1) *= *new_c2;
////        }
////        if(c1.is_number()) {
////            if (c1._is_transposed && !c1.is_matrix()) {//TODO check when matrix
////                auto new_c2 = copy(c2);
////                new_c2->_is_vector = true;
////                auto res = func_(c1);
////                res._dim[1] = c2._dim[0];
////                res *= move(*new_c2);
////                res._is_vector = false;
////                res._is_transposed = false;
////                delete new_c2;
////                return res;
////            }
////            else {
////                return func_(c1) *= c2;
////            }
////        }
////        else {
//            //            if (c1._is_transposed) {
//            ////                auto new_c2 = copy(c2);
//            ////                new_c2->_is_vector = true;
//            ////                auto res = func_(c1) *= move(*new_c2);
//            ////                auto res = ;
//            ////                delete new_c2;
//            //                return func_(c1) *= c2;
//            //            }
//            //            if (c2.is_var()) {
//            //                return func_(c2) *= c1;
//            //            }
//            //            if(c1.is_function()){
//            //                return (*(func_*)&c1) *= c2;
//            //            }
//            //            else {
//            return func_(c1) *= c2;
//            //            }
//        }
//    }
//
//    //    func_ operator*(func_&& f, const constant_& c){
//    //        return f *= c;
//    //    }
//    //
//    //    func_ operator*(const constant_& c, func_&& f){
//    //        return f *= c;
//    //    }
//
//    func_ operator/(const constant_& c1, const constant_& c2){
//        return func_(c1) /= c2;
//    }
//
//    //    func_ operator/(func_&& f, const constant_& c){
//    //        return f /= c;
//    //    }
//
//    //
//    //func_ operator+(const func_& f1, const func_& f2){
//    //    return func_(f1) += f2;
//    //}
//    //
//    //func_ operator+(func_&& f1, const func_& f2){
//    //    return f1 += f2;
//    //}
//    //
//    //func_ operator+(const func_& f1, func_&& f2){
//    //    return f2 += f1;
//    //}
//    //
//    //
//    //func_ operator-(const func_& f1, const func_& f2){
//    //    return func_(f1) -= f2;
//    //}
//    //
//    //func_ operator-(func_&& f1, const func_& f2){
//    //    return f1 -= f2;
//    //}
//    //
//    //func_ operator-(const func_& f1, func_&& f2){
//    //    return (f2 *= -1) += f1;
//    //}
//    //
//    constant_* add(constant_* c1, const func_& f){
//        switch (c1->get_type()) {
//            case binary_c: {
//                auto res = new func_(f);
//                *res += ((constant<bool>*)c1)->eval();
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case short_c: {
//                auto res = new func_(f);
//                *res += ((constant<short>*)c1)->eval();
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case integer_c: {
//                auto res = new func_(f);
//                *res += ((constant<int>*)c1)->eval();
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case float_c: {
//                auto res = new func_(f);
//                *res += ((constant<float>*)c1)->eval();
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case double_c: {
//                auto res = new func_(f);
//                *res += ((constant<double>*)c1)->eval();
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case long_c: {
//                auto res = new func_(f);
//                *res += ((constant<long double>*)c1)->eval();
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case complex_c:{
//                auto res = new func_(f);
//                *res += *((constant<Cpx>*)c1);
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case par_c:{
//                auto res = new func_(f);
//                if (f.is_constant()) {
//                    res->insert(true, constant<double>(1), *(param_*)c1);
//                }
//                else{
//                    auto cst = res->get_cst();
//                    cst = add(cst, *c1);
//                }
//                delete c1;
//                return c1 = res;
//                break;
//            }
//            case var_c:{
//                auto res = new func_(f);
//                delete c1;
//                res->insert(true, constant<double>(1), *(param_*)c1);
//                return c1 = res;
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(plus_, make_shared<func_>((*c1)), make_shared<func_>(f));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(plus_, make_shared<func_>((*c1)), make_shared<func_>(f));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//                auto res = new func_(f);
//                *res += *((func_*)c1);
//                delete c1;
//                return c1 = res;
//                break;        }
//
//            default:
//                break;
//        }
//        return c1;
//    }
//
//    constant_* add(constant_* c1, const constant_& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c2.get_type()) {
//            case binary_c: {
//                return add(c1, *(constant<bool>*)&c2);
//                break;
//            }
//            case short_c: {
//                return add(c1, *(constant<short>*)&c2);
//                break;
//            }
//            case integer_c: {
//                return add(c1, *(constant<int>*)&c2);
//                break;
//            }
//            case float_c: {
//                return add(c1, *(constant<float>*)&c2);
//                break;
//            }
//            case double_c: {
//                return add(c1, *(constant<double>*)&c2);
//                break;
//            }
//            case long_c: {
//                return add(c1, *(constant<long double>*)&c2);
//                break;
//            }
//            case complex_c: {
//                return add(c1, *(constant<Cpx>*)&c2);
//                break;
//            }
//            case par_c:{
//                auto pc2 = (param_*)(&c2);
//                switch (pc2->get_intype()) {
//                    case binary_:
//                        return add(c1, *(param<bool>*)pc2);
//                        break;
//                    case short_:
//                        return add(c1, *(param<short>*)pc2);
//                        break;
//                    case integer_:
//                        return add(c1, *(param<int>*)pc2);
//                        break;
//                    case float_:
//                        return add(c1, *(param<float>*)pc2);
//                        break;
//                    case double_:
//                        return add(c1, *(param<double>*)pc2);
//                        break;
//                    case long_:
//                        return add(c1, *(param<long double>*)pc2);
//                        break;
//                    case complex_:
//                        return add(c1, *(param<Cpx>*)pc2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case var_c:{
//                auto pc2 = (param_*)(&c2);
//                switch (pc2->get_intype()) {
//                    case binary_:
//                        return add(c1, *(var<bool>*)pc2);
//                        break;
//                    case short_:
//                        return add(c1, *(var<short>*)pc2);
//                        break;
//                    case integer_:
//                        return add(c1, *(var<int>*)pc2);
//                        break;
//                    case float_:
//                        return add(c1, *(var<float>*)pc2);
//                        break;
//                    case double_:
//                        return add(c1, *(var<double>*)pc2);
//                        break;
//                    case long_:
//                        return add(c1, *(var<long double>*)pc2);
//                        break;
//                    case complex_:
//                        return add(c1, *(var<Cpx>*)pc2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(plus_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(plus_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//                return add(c1, *(func_*)&c2);
//                break;
//            }
//            default:
//                break;
//        }
//        return nullptr;
//    }
//
//    constant_* substract(constant_* c1, const constant_& c2){ /**< substracts c2 from c1, updates its type and returns the result **/
//        switch (c2.get_type()) {
//            case binary_c: {
//                return substract(c1, *(constant<bool>*)&c2);
//                break;
//            }
//            case short_c: {
//                return substract(c1, *(constant<short>*)&c2);
//                break;
//            }
//            case integer_c: {
//                return substract(c1, *(constant<int>*)&c2);
//                break;
//            }
//            case float_c: {
//                return substract(c1, *(constant<float>*)&c2);
//                break;
//            }
//            case double_c: {
//                return substract(c1, *(constant<double>*)&c2);
//                break;
//            }
//            case long_c: {
//                return substract(c1, *(constant<long double>*)&c2);
//                break;
//            }
//            case complex_c: {
//                return substract(c1, *(constant<Cpx>*)&c2);
//                break;
//            }
//            case par_c:{
//                auto pc2 = (param_*)(&c2);
//                switch (pc2->get_intype()) {
//                    case binary_:
//                        return substract(c1, *(param<bool>*)pc2);
//                        break;
//                    case short_:
//                        return substract(c1, *(param<short>*)pc2);
//                        break;
//                    case integer_:
//                        return substract(c1, *(param<int>*)pc2);
//                        break;
//                    case float_:
//                        return substract(c1, *(param<float>*)pc2);
//                        break;
//                    case double_:
//                        return substract(c1, *(param<double>*)pc2);
//                        break;
//                    case long_:
//                        return substract(c1, *(param<long double>*)pc2);
//                        break;
//                    case complex_:
//                        return substract(c1, *(param<Cpx>*)pc2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case var_c:{
//                auto pc2 = (param_*)(&c2);
//                switch (pc2->get_intype()) {
//                    case binary_:
//                        return substract(c1, *(var<bool>*)pc2);
//                        break;
//                    case short_:
//                        return substract(c1, *(var<short>*)pc2);
//                        break;
//                    case integer_:
//                        return substract(c1, *(var<int>*)pc2);
//                        break;
//                    case float_:
//                        return substract(c1, *(var<float>*)pc2);
//                        break;
//                    case double_:
//                        return substract(c1, *(var<double>*)pc2);
//                        break;
//                    case long_:
//                        return substract(c1, *(var<long double>*)pc2);
//                        break;
//                    case complex_:
//                        return substract(c1, *(var<Cpx>*)pc2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(minus_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(minus_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case func_c: {
//                auto f2 = (func_*)&c2;
//                if (f2->is_number()) {
//                    return substract(c1, *f2->get_cst());
//                }
//                auto f = new func_(*c1);
//                delete c1;
//                *f -= *f2;
//                return c1 = (constant_*)f;
//                break;
//                break;
//            }
//        }
//        return nullptr;
//    }
//
//    constant_* multiply(constant_* c1, const constant_& c2){ /**< adds c2 to c1, updates its type and returns the result **/
//        switch (c2.get_type()) {
//            case binary_c: {
//                return multiply(c1, *(constant<bool>*)&c2);
//                break;
//            }
//            case short_c: {
//                return multiply(c1, *(constant<short>*)&c2);
//                break;
//            }
//            case integer_c: {
//                return multiply(c1, *(constant<int>*)&c2);
//                break;
//            }
//            case float_c: {
//                return multiply(c1, *(constant<float>*)&c2);
//                break;
//            }
//            case double_c: {
//                return multiply(c1, *(constant<double>*)&c2);
//                break;
//            }
//            case long_c: {
//                return multiply(c1, *(constant<long double>*)&c2);
//                break;
//            }
//            case complex_c: {
//                return multiply(c1, *(constant<Cpx>*)&c2);
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case par_c:{
//                auto pc2 = (param_*)(&c2);
//                switch (pc2->get_intype()) {
//                    case binary_:
//                        return multiply(c1, *(param<bool>*)pc2);
//                        break;
//                    case short_:
//                        return multiply(c1, *(param<short>*)pc2);
//                        break;
//                    case integer_:
//                        return multiply(c1, *(param<int>*)pc2);
//                        break;
//                    case float_:
//                        return multiply(c1, *(param<float>*)pc2);
//                        break;
//                    case double_:
//                        return multiply(c1, *(param<double>*)pc2);
//                        break;
//                    case long_:
//                        return multiply(c1, *(param<long double>*)pc2);
//                        break;
//                    case complex_:
//                        return multiply(c1, *(param<Cpx>*)pc2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case var_c:{
//                auto pc2 = (param_*)(&c2);
//                switch (pc2->get_intype()) {
//                    case binary_:
//                        return multiply(c1, *(var<bool>*)pc2);
//                        break;
//                    case short_:
//                        return multiply(c1, *(var<short>*)pc2);
//                        break;
//                    case integer_:
//                        return multiply(c1, *(var<int>*)pc2);
//                        break;
//                    case float_:
//                        return multiply(c1, *(var<float>*)pc2);
//                        break;
//                    case double_:
//                        return multiply(c1, *(var<double>*)pc2);
//                        break;
//                    case long_:
//                        return multiply(c1, *(var<long double>*)pc2);
//                        break;
//                    case complex_:
//                        return multiply(c1, *(var<Cpx>*)pc2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case func_c: {
//                auto f = new func_(c2);
//                *f *= *c1;
//                delete c1;
//                return c1 = (constant_*)f;
//                break;
//            }
//            default:
//                break;
//        }
//        return nullptr;
//    }
//
//
//    constant_* divide(constant_* c1, const constant_& c2){
//        switch (c2.get_type()) {
//            case binary_c: {
//                return divide(c1, *(constant<bool>*)&c2);
//                break;
//            }
//            case short_c: {
//                return divide(c1, *(constant<short>*)&c2);
//                break;
//            }
//            case integer_c: {
//                return divide(c1, *(constant<int>*)&c2);
//                break;
//            }
//            case float_c: {
//                return divide(c1, *(constant<float>*)&c2);
//                break;
//            }
//            case double_c: {
//                return divide(c1, *(constant<double>*)&c2);
//                break;
//            }
//            case long_c: {
//                return divide(c1, *(constant<long double>*)&c2);
//                break;
//            }
//            case complex_c: {
//                return divide(c1, *(constant<Cpx>*)&c2);
//                break;
//            }
//            case uexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            case bexp_c: {
//                auto res = new bexpr(product_, make_shared<func_>((*c1)), make_shared<func_>((c2)));
//                delete c1;
//                c1 = (constant_*)res;
//                return c1;
//                break;
//            }
//            default:{
//                auto f = new func_(*c1);
//                *f /= func_(c2);
//                delete c1;
//                return c1 = (constant_*)f;
//                break;
//            }
//        }
//        return nullptr;
//    }
//
//
//    shared_ptr<func_> func_::get_stored_derivative(const string& vid) const{
//        auto it = _dfdx->find(vid);
//        if (it!=_dfdx->end()) {
//            return it->second;
//        }
//        else {
//            throw invalid_argument("No derivatives stored!\n");
//        }
//    }
//
//    func_ func_::get_dfdx(const param_& v){
//        _dfdx->clear();
//        compute_derivatives();
//        return *get_stored_derivative(v._name).get();
//    }
//
//    func_* func_::compute_derivative(const param_ &v){
//        auto vid = v._name;
//        if(_dfdx->count(vid)!=0){
//            return _dfdx->at(vid).get();
//        }
//
//        auto df = new func_(get_derivative(v));
//        df->allocate_mem();
//        (*_dfdx)[vid] = shared_ptr<func_>(df);
//        DebugOff( "First derivative with respect to " << v.get_name(false,false) << " = " << df->to_str() << endl);
//        return df;
//    }
//
//    void func_::compute_derivatives(){ /**< Computes and stores the derivative of f with respect to all variables. */
//        size_t vid = 0, vjd = 0;
//        param_* vi;
//        param_* vj;
//        DebugOff( "Computing derivatives for " << to_str() << endl);
//        for (auto &vp: *_vars) {
//            vi = vp.second.first.get();
//            vid = vi->get_id();
//            auto vi_name = vp.first;
//            auto df = compute_derivative(*vi);
//            //            if (is_nonlinear()) {
//            DebugOff( "First derivative with respect to " << vp.first << " = " << df->to_str() << endl);
//            //            }
//            for (auto &vp2: *df->_vars) {
//                vj = vp2.second.first.get();
//                vjd = vj->get_id();
//                auto vj_name = vp2.first;
//                if (vi_name.compare(vj_name) <= 0) { //only store lower left part of hessian matrix since it is symmetric.
//                    auto d2f = df->compute_derivative(*vj);
//                    DebugOff( "Second derivative with respect to " << vp2.first << " and " << vp.first << " = " << d2f->to_str() << endl);
//                    //                d2f->print();
//                }
//            }
//
//        }
//    };
//
    /* Build Compositions of positive integer n given k positive integers. All credits go to Dr. Martin von Gagern */
    vector<vector<int>> build_compositions(int k, int n) {
        int pos = 0, last = n - 1;
        vector<vector<int>> res;
        vector<int> a;
        a.resize(n);
        a[0] = k;
        res.push_back(a);
        while (true) {
            if (pos != last) {
                --a[pos];
                ++pos;
                a[pos] = 1;
            }
            else {
                if (a[last] == k)
                    return res;
                for (--pos; a[pos] == 0; --pos);
                --a[pos];
                int tmp = 1 + a[last];
                ++pos;
                a[last] = 0;
                a[pos] = tmp;
            }
            res.push_back(a);
        }
    }

    /* Build all monomials of degree d in dimension n (size of vars) */
    vector<pterm> func_::get_monomials(unsigned d){
        vector<pterm> res;
        int n = _vars->size();
        auto comp = build_compositions(d, n);
        auto m = comp.size();
        res.resize(m);
        for (auto i = 0; i<m; i++) {
            auto l = make_shared<list<pair<shared_ptr<param_>, int>>>();
            auto row = comp[i];
            auto v_it = _vars->begin();
            for (auto j = 0; j<n; j++) {
                auto exp = row[j];
                auto v = v_it->second.first;
                if (exp>0) {
                    l->push_back(make_pair<>(v, exp));
                }
                v_it++;
            }
            res[i] = pterm(true, make_shared<constant<>>(1), l);
            cout << res[i].to_str(0, 1) << " | ";
        }
        cout << endl;
        return res;
    }
//
//
//
//    void func_::reset_val(){
//        _evaluated = false;
//        if (_expr) {
//            _expr->reset_val();
//        }
//    }
//
//    void func_::reindex(param_* v){
//
//    }
//
//    void func_::replace(param_* v, func_& f){
//        auto vin = get_var(*v->_vec_id);
//        f.reindex(v);
//        //index f here!
//        for (auto &l_p:*_lterms) {
//            auto &lt = l_p.second;
//            if (lt._p==vin) {
//                //replace here.
//            }
//        }
//
//        if (_expr) {
//            if (_expr->is_uexpr()) {
//                auto ue = (uexpr*)_expr.get();
//                ue->_son->replace(v,f);
//            }
//            else {
//                auto be = (bexpr*)_expr.get();
//                be->_lson->replace(v,f);
//                be->_rson->replace(v,f);
//            }
//        }
//        _vars->erase(v->_name);
//        DebugOn(to_str());
//        //        print();
//    }
//
    shared_ptr<param_> func_::get_var(size_t vid) const{
        for (auto &v_p:*_vars) {
            if (*v_p.second.first->_vec_id == vid) {
                return v_p.second.first;
            }
        }
        return nullptr;
    }
//
    bool func_::has_var(const param_& v) const{
        return get_var(v.get_name(true, false))!=nullptr;
    }
        
//
//
//    bool func_::has_var(const string& name) const{
//        return _vars->count(name)>0;
//    }
//
//    func_ func_::get_derivative(const param_ &v) const{
//        func_ res;
//        if(!has_var(v)){
//            return res;
//        }
//        for (auto &lt: *_lterms) {
//            if (*lt.second._p == v) {
//                auto coef = copy(*lt.second._coef);
//                if ((coef->_is_vector && coef->_is_transposed)) {
//                    coef->transpose();
//                }
//                res = move(*coef);
//                delete coef;
//                if(!lt.second._sign){
//                    res *= -1;
//                }
//            }
//        }
//        for (auto &lt: *_qterms) {
//            if (*lt.second._p->first == v) {
//                auto coef = copy(*lt.second._coef);
//                if (coef->_is_vector && coef->_is_transposed) {
//                    coef->transpose();
//                }
//                if(lt.second._sign) {
////                    if(lt.second._c_p1_transposed){
////                        res += (*coef*(*lt.second._p->second)).tr();
////                    }
////                    else {
//                        res += *coef*(*lt.second._p->second);
////                    }
//                }
//                else {
////                    if(lt.second._c_p1_transposed){
////                        res -= (*coef*(*lt.second._p->second)).tr();
////                    }
////                    else {
//                        res -= *coef*(*lt.second._p->second);
////                    }
//                }
//                delete coef;
//                if (lt.second._coef->_is_vector || lt.second._coef->is_matrix()) {
//                    res._is_vector = true;
//                    //                    res._nb_instances = lt.second._coef->get_dim();
//                }
//            }
//            if (*lt.second._p->second == v) {
//                auto coef = copy(*lt.second._coef);
//                if (coef->_is_vector && coef->_is_transposed) {
//                    coef->transpose();
//                }
//                if(lt.second._sign) {
//                    res += *coef*(*lt.second._p->first);
//                }
//                else {
//                    res -= *coef*(*lt.second._p->first);
//                }
//                delete coef;
//                if (lt.second._coef->_is_vector || lt.second._coef->is_matrix()) {
//                    res._is_vector = true;
//                }
//            }
//
//            //CHECK THIS
//            //            if (lt.second._coef->_is_vector || lt.second._coef->_is_matrix) {
//            //                res._is_vector = true;
//            //                res._nb_instances = lt.second._coef->get_dim();
//            //                res._dim[0] = lt.second._coef->_dim[0]; res._dim[1] = lt.second._coef->_dim[1];
//            //            }
//        }
//        for (auto &lt: *_pterms) {
//            for (auto &p: *lt.second._l) {
//                if (*p.first == v) {
//                    func_ pterm = constant<double>(1);
//                    if (!lt.second._sign) {
//                        pterm = constant<double>(-1);
//                    }
//                    auto expo = p.second;
//                    if (expo > 1) {
//                        pterm *= expo;
//                        pterm *= *lt.second._coef;
//                        pterm *= (*p.first);
//                        for (int i = 1; i<expo-1; i++) {
//                            pterm *= *p.first;
//                        }
//                    }
//                    for (auto &p2: *lt.second._l) {
//                        if (p2!=p) {
//                            func_ pterm2(*p2.first);
//                            for (int i = 1; i<p2.second; i++) {
//                                pterm2 *= *p2.first;
//                            }
//                            pterm *= pterm2;
//                        }
//                    }
//                    res += *lt.second._coef*pterm;
//                }
//            }
//        }
//        if (!_expr) {
//            //            res.untranspose();
//            return res;
//        }
//        else { // f is a composition of functions
//            res += _expr->get_derivative(v);
//            //            res.untranspose();
//            return res;
//        }
//    }
//
//    double func_::get_val(size_t i) const{
//        //        if (is_constant()) {
//        //            if(eval(i)!=_val->at(i)){
//        //                throw invalid_argument("error");
//        //            }
//        //        }
//        //        if (is_number()) {
//        //            return _val->at(0);
//        //        }
//        //        else {
//        //            if (i>=_val->size()) {
//        //                throw invalid_argument("error");
//        //            }
//        if (_val->size()<=i){
//            throw invalid_argument("Func get_val out of range");
//        }
//        return _val->at(i);
//        //        }
//    }
//
//    double func_::get_val(size_t i, size_t j) const{
//
//        //        if () {
//        //            return _val->at(0);
//        //        }
//        //        if (!_is_matrix) {
//        //            //            if (_is_transposed) {
//        //            //                return get_val(i);
//        //            //            }
//        //            throw invalid_argument("get_val i,j");
//        //            return get_val(j);
//        //        }
//        //        if (is_number()) {
//        //            return _val->at(0);
//        //        }
//        //        else {
//        if (_is_transposed) {
//            //                if (j*_dim[0]+i>=_val->size()) {
//            //                    throw invalid_argument("error");
//            //                }
//            return _val->at(j*_dim[0]+i);
//        }
//        //            if (i*_dim[1]+j>=_val->size()) {
//        //                throw invalid_argument("error");
//        //            }
//        return _val->at(i*_dim[1]+j);
//        //        }
//    }
//
//    //    double func_::force_eval(size_t i){
//    //
//    //        double res = 0;
//    //        for (auto &pair:*_pterms) {
//    //            res += pair.second.eval(i);
//    //        }
//    //        for (auto &pair:*_qterms) {
//    //            res += pair.second.eval(i);
//    //        }
//    //        for (auto &pair:*_lterms) {
//    //            res += pair.second.eval(i);
//    //        }
//    //        res += eval(_cst,i);
//    //        if (_expr) {
//    //            res += _expr->eval(i);
//    //        }
//    //        return res;
//    //    }
//
//    double func_::eval(size_t i, size_t j){
//        //        if (_val->size()<_nb_instances) {
//        //            _val->resize(_nb_instances);
//        //        }
//        if (!is_matrix() && (!_indices || !_indices->_ids || _indices->_ids->size()==1)) {
//            //            if (!is_matrix()) {
//            //            if (_is_transposed) {
//            //                return eval(i);
//            //            }
//            return eval(j);
//        }
//
//        if (is_constant() && _evaluated) {
//            if (is_number()) {
//                return _val->at(0);
//            }
//            else {
//                if (_is_transposed) {
//                    return _val->at(j*_dim[0]+i);
//                }
//                return _val->at(i*_dim[1]+j);
//            }
//        }
//        double res = 0;
//        for (auto &pair:*_pterms) {
//            res += pair.second.eval(i,j);
//        }
//        for (auto &pair:*_qterms) {
//            res += pair.second.eval(i,j);
//        }
//        for (auto &pair:*_lterms) {
//            res += pair.second.eval(i,j);
//        }
//        res += t_eval(_cst,i,j);
//        if (_expr) {
//            res += _expr->eval(i,j);
//        }
//        if (is_number()) {
//            _val->at(0) = res;
//        }
//        else {
//            if (is_constant()) {
//                if (_is_transposed) {
//                    if(_dim[0]*j + i ==_val->size()-1){
//                        _evaluated = true;
//                    }
//
//                }
//                else {
//                    if(_dim[1]*i + j ==_val->size()-1){
//                        _evaluated = true;
//                    }
//                }
//
//            }
//            set_val(i,j,res);
//        }
//        return res;
//    }
//
//    double func_::eval(size_t i){
//                if (_val->size()<=i) {
//                    allocate_mem();
//                }
//        //        if (!_val) {
//        //            throw invalid_argument("_val not defined for function.\n");
//        //        }
//        if (is_constant() && _evaluated) {
//            if (is_number()) {
//                return _val->at(0);
//            }
//            //            if (i>=_val->size()) {
//            //                throw invalid_argument("error");
//            //            }
//            //            if (_val->at(i) != force_eval(i)) {
//            //                throw invalid_argument("error");
//            //            }
//            if (_val->size()<=i){
//                throw invalid_argument("Func eval out of range");
//            }
//            return _val->at(i);
//        }
//        double res = 0;
//        for (auto &pair:*_pterms) {
//            res += pair.second.eval(i);
//        }
//        for (auto &pair:*_qterms) {
//            res += pair.second.eval(i);
//        }
//        for (auto &pair:*_lterms) {
//            res += pair.second.eval(i);
//        }
//        res += t_eval(_cst,i);
//        if (_expr) {
//            //            if (_expr->is_uexpr()) {
//            //                auto ue = (uexpr*)_expr.get();
//            ////                if (ue->_son->is_constant()) {
//            ////                _nb_instances = max(_nb_instances, ue->_son->_nb_instances);
//            ////                    _val->resize(max(_val->size(),ue->_son->_nb_instances));//TODO is this necessary?
//            ////                }
//            //
//            //            }
//            //            else {
//            //                auto be = (bexpr*)_expr.get();
//            //                if(!be->is_inner_product()) {
//            //                    _nb_instances = max(_nb_instances, max(be->_lson->_nb_instances,be->_rson->_nb_instances));
//            //                }
//            ////                if (be->_lson->is_constant() && be->_rson->is_constant()) {
//            ////                    _val->resize(max(_val->size(),max(be->_lson->_nb_instances,be->_rson->_nb_instances)));
//            ////                }
//            //
//            //            }
//            res += _expr->eval(i);
//        }
//        if (is_number()) {
//            _val->at(0) = res;
//            _evaluated = true;//TODO fix this
//        }
//        else {
//            //            if (i>=_val->size()) {
//            //                throw invalid_argument("error");
//            //            }
//            if (is_constant() && i==get_dim()-1) {
//                _evaluated = true;
//            }
//            if (_val->size()<=i){
//                throw invalid_argument("Param eval out of range");
//            }
//            _val->at(i) = res;
//        }
//        return res;
//    }
//
//    string func_::to_str() const{ return string();};
//        string str;
//        int ind = 0;
//        string sign = " + ";
//        //        for (int inst = 0; inst < _nb_instances ; inst++) {
//        ind = 0;
//
//        for (auto &pair:*_pterms) {
//            str += pair.second.to_str(ind++);
//        }
//        if (!_pterms->empty() && (!_qterms->empty() || !_lterms->empty())) {
//            str += " + ";
//        }
//        ind = 0;
//        for (auto &pair:*_qterms) {
//            str += pair.second.to_str(ind++);
//        }
//        if (!_qterms->empty() && !_lterms->empty()) {
//            str += " + ";
//        }
//        ind = 0;
//        for (auto &pair:*_lterms) {
//            str += pair.second.to_str(ind++);
//        }
//        if (_cst->is_number()) {
//            auto val = poly_to_str(_cst);
//            if (val.front()=='-') {
//                str += " - " + val.substr(1);
//            }
//            else if (val != "0"){
//                if (!_pterms->empty() || !_qterms->empty() || !_lterms->empty()) {
//                    str += " + ";
//                }
//                str += val;
//            }
//        }
//        else {
//            if (!_pterms->empty() || !_qterms->empty() || !_lterms->empty()) {
//                str += " + ";
//            }
//            str += "(";
//            str += poly_to_str(_cst);
//            str += ")";
//        }
//        if (_expr && (!_pterms->empty() || !_qterms->empty() || !_lterms->empty() || !_cst->is_zero())) {
//            str += " + ";
//        }
//        if (_expr) {
//            str += _expr->get_str();
//        }
//        if (_is_vector && (is_number() || _vars->size()>1 || _params->size()>1)) {
//            str = "[" + str +"]";
//        }
//        if (_is_transposed) {
//            str += "\u1D40";
//        }
//        if(str.size()==0){
//            str = "0";
//        }
//        return str;
//    }
//
//
    void func_::print_symbolic(bool endline, bool display_input){
        string str;
        if (display_input) {
            if (is_constant()) {
                cout << " (Constant";
            }
            else if (is_linear()) {
                cout << " (Linear";
            }
            else if (is_convex()) {
                cout << " (Convex";
            }
            else if (is_concave()){
                cout << " (Concave";
            }
            else {
                cout << " (Unknown";
            }
            if (is_complex()) {
                cout << " Complex) : ";
            }
            else {
                cout << ") : ";
            }
            if (!_embedded && !is_constant()) {
                str += "f(";
                for (auto pair_it = _vars->begin(); pair_it != _vars->end();) {
                    str += pair_it->second.first->get_name(false,false);
                    if (next(pair_it) != _vars->end()) {
                        str += ",";
                    }
                    pair_it++;
                }
                str += ") = ";
            }
        }
        str += to_str();
        _to_str = str;
        cout << this->_to_str;
        if (endline)
            cout << endl;
    }

    

    
    
//
    
//
//    void func_::print(size_t index) {
//        cout << to_str(index);
//    }
//
//    void func_::print(){
//        string str;
//        if (is_constant()) {
//            str += " (Constant) : ";
//        }
//        else if (is_complex()) {
//            cout << " (Complex) : ";
//        }
//        else if (is_linear()) {
//            str += " (Linear) : ";
//        }
//        else if (is_convex()) {
//            str += " (Convex) : ";
//        }
//        else if (is_concave()){
//            str += " (Concave) : ";
//        }
//        else {
//            str += " (Unknown) : ";
//        }
//        if (!_embedded && !is_constant()) {
//            str += "f(";
//            for (auto pair_it = _vars->begin(); pair_it != _vars->end();) {
//                //                    if (!pair_it->second.first->_is_vector) {
//                str += pair_it->second.first->get_name(false,true);
//                //                    str += pair_it->second.first->get_name(false,false)+"[";
//                //                    str += to_string(pair_it->second.first->get_id_inst())+"]";
//                if (next(pair_it) != _vars->end()) {
//                    str += ",";
//                }
//                //                }
//                pair_it++;
//            }
//            str += ") = ";
//        }
//        auto space_size = str.size();
//        auto nb_inst = _dim[0];
//        allocate_mem();
//        for (size_t inst = 0; inst<nb_inst; inst++) {
//            eval(inst);
//            if (inst>0) {
//                str.insert(str.end(), space_size, ' ');
//            }
//            str += to_str(inst);
//            str += "\n";
//        }
//        str += "\n";
//        cout << str;
//    }
//
//
//    void func_::update_to_str(bool input){
//        _to_str = to_str();
//    }


    /**
     Returns a pointer to the constant part of the function.
     @return a pointer to the constant part of the function.
     */
    shared_ptr<constant_> func_::get_cst() const{
        return _cst;
    }

    /**
     Returns a pointer to the variable matching the name provided.
     @param[in] name variable name.
     @return a pointer to the variable matching the name provided.
     */
    shared_ptr<param_> func_::get_var(const string& name) const{
        if (_vars->empty()) {
            return nullptr;
        }
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }

    /**
     Returns a pointer to the parameter matching the name provided.
     @param[in] name variable name.
     @return a pointer to the parameter matching the name provided.
     */
    shared_ptr<param_> func_::get_param(string name) const{
        auto pair_it = _params->find(name);
        if (pair_it==_params->end()) {
            return nullptr;
        }
        else {
            return get<1>(*pair_it).first;
        }
    }

    void func_::add_var(shared_ptr<param_> v, int nb){/**< Inserts the variable in this function input list. nb represents the number of occurences v has. WARNING: Assumes that v has not been added previousely!*/
        if (_vars->count(v->get_name(false,false))!=0) {
            throw invalid_argument("In function add_var(v,nb): Variable already contained in function");
        }
        _vars->insert(make_pair<>(v->get_name(false,false), make_pair<>(v, nb)));
        if(v->is_double_indexed()){
            _indices = v->_indices;
        }
        if (v->_is_vector) {// i.e., it appears in a sum
            if (v->is_matrix()) {
                if (v->_is_transposed) {
                    _nb_vars += v->get_dim(0);
                }
                _nb_vars += v->get_dim(1);
            }
            else {
                _nb_vars += v->get_dim();
            }
        }
        else {
            _nb_vars++;
        }
    }
    
    func<double> get_mag(constant_* c){
        if(c->is_double()){
            func<> res = *((constant<double>*)c);
            return res;
        }
        if(c->is_complex()){
            func<> res = std::abs(((constant<Cpx>*)c)->eval());
            if(c->_is_transposed){
                res.transpose();
            }
            return res;
        }
        if(c->is_param() || c->is_var()){
            auto p = (param_*)c;
            if(p->_is_imag || p->_is_angle || p->_is_real || p->_is_sqrmag){
                throw invalid_argument("unsupported");
            }
            if(p->get_intype()==complex_){
                if(p->_mag){
                    if(p->is_var()){
                        return *((var<>*)p->_mag.get());
                    }
                    else {
                        return *((param<>*)p->_mag.get());
                    }
                }
                else {
                    return 0;
                }
            }
            else if(p->get_intype()==double_){
                if(p->is_var()){
                    return *((var<>*)p);
                }
                else {
                    return *((param<>*)p);
                }
            }
        }
        if(c->is_function()){
            auto f = (func_*)c;
            if(f->get_return_type()==double_){
                return *(func<>*)f;
            }
            if(f->get_return_type()==complex_){
                return get_mag_ang(*(func<Cpx>*)f).first;
            }
        }
        throw invalid_argument("unsupported");
    }
    
    func<double> get_real(constant_* c){
        if(c->is_double()){
            func<> res = ((constant<double>*)c)->eval();
            if(c->_is_transposed){
                res.transpose();
            }
            return res;
        }
        if(c->is_complex()){
            func<> res = ((constant<Cpx>*)c)->eval().real();
            if(c->_is_transposed){
                res.transpose();
            }
            return res;
        }
        if(c->is_param() || c->is_var()){
            auto p = (param_*)c;
            if(p->_is_imag || p->_is_angle || p->_is_real || p->_is_sqrmag){
                throw invalid_argument("unsupported");
            }
            if(p->get_intype()==complex_){
                if(p->_real){
                    if(p->is_var()){
                        return *((var<>*)p->_real.get());
                    }
                    else {
                        return *((param<>*)p->_real.get());
                    }
                }
                else {
                    return 0;
                }
            }
            else if(p->get_intype()==double_){
                if(p->is_var()){
                    return *((var<>*)p);
                }
                else {
                    return *((param<>*)p);
                }
            }
        }
        if(c->is_function()){
            auto f = (func_*)c;
            if(f->get_return_type()==double_){
                return *(func<>*)f;
            }
            if(f->get_return_type()==complex_){
                return get_real_imag(*(func<Cpx>*)f).first;
            }
        }
        throw invalid_argument("unsupported");
    }
    
    func<double> get_ang(constant_* c){
        if(c->is_double()){
            return std::atan2(0,((constant<double>*)c)->eval());
        }
        if(c->is_complex()){
            return arg(((constant<Cpx>*)c)->eval());
        }
        if(c->is_param() || c->is_var()){
            auto p = (param_*)c;
            if(p->_is_angle){
                if(p->is_var()){
                    return *((var<>*)p);
                }
                else {
                    return *((param<>*)p);
                }
            }
            if(p->_is_imag || p->_is_real || p->_is_sqrmag){
                throw invalid_argument("unsupported");
            }
            if(p->get_intype()==complex_){
                if(p->_ang){
                    if(p->is_var()){
                        if(p->_is_conjugate){
                            func<> res = (*((var<>*)p->_ang.get()));
                            res.reverse_sign();
                            return res;
                        }
                        else {
                            return *((var<>*)p->_ang.get());
                        }
                    }
                    else {
                        if(p->_is_conjugate){
                            func<> res = (*((param<>*)p->_ang.get()));
                            res.reverse_sign();
                            return res;
                        }
                        else {
                            return *((param<>*)p->_ang.get());
                        }
                    }
                }
                else {
                    return 0;
                }
            }
            else {
                return 0;
            }
        }
        if(c->is_function()){
            auto f = (func_*)c;
            if(f->get_return_type()==complex_){
                return get_mag_ang(*(func<Cpx>*)f).second;
            }
            return 0;
        }
        throw invalid_argument("unsupported");
    }
    
    func<double> get_imag(constant_* c){
        if(c->is_double()){
            return 0;
        }
        if(c->is_complex()){
            return ((constant<Cpx>*)c)->eval().imag();
        }
        if(c->is_param() || c->is_var()){
            auto p = (param_*)c;
            if(p->_is_imag || p->_is_angle || p->_is_real || p->_is_sqrmag){
                throw invalid_argument("unsupported");
            }
            if(p->get_intype()==complex_){
                if(p->_imag){
                    if(p->is_var()){
                        if(p->_is_conjugate){
                            func<> res = (*((var<>*)p->_imag.get()));
                            res.reverse_sign();
                            return res;
                        }
                        else {
                            return *((var<>*)p->_imag.get());
                        }
                    }
                    else {
                        if(p->_is_conjugate){
                            func<> res = (*((param<>*)p->_imag.get()));
                            res.reverse_sign();
                            return res;
                        }
                        else {
                            return *((param<>*)p->_imag.get());
                        }
                    }
                }
                else {
                    return 0;
                }
            }
            else {
                return 0;
            }
        }
        if(c->is_function()){
            auto f = (func_*)c;
            if(f->get_return_type()==complex_){
                return get_real_imag(*(func<Cpx>*)f).second;
            }
            return 0;
        }
        throw invalid_argument("unsupported");
    }
    
    
    void func_::add_param(shared_ptr<param_> p, int nb){/**< Inserts the parameter in this function input list. WARNING: Assumes that p has not been added previousely!*/
        if (_params->count(p->get_name(false,false))!=0) {
            throw invalid_argument("In function add_param(v,nb): parameter already contained in function");
        }
//        set_max_dim(*p);
        _params->insert(make_pair<>(p->get_name(false,false), make_pair<>(p, nb)));
    }
    
//    template<>
//    double func<double>::eval_double(size_t i) {
////        if(is_zero()){
////            return _range->first;
////        }
//        if (_vars->empty() && _evaluated) {
//            if (_params->empty()){
//                return _val->at(0);
//            }
//            return _val->at(i);
//        }
//        double res = 0;
////        if(!_cst->is_zero())
//            res += eval_cst(i);
////        if(!_lterms->empty()){
//            for (auto &pair:*_lterms) {
//                if ((pair.second._coef->_is_transposed || pair.second._coef->is_matrix() || (pair.second._p->is_indexed() && pair.second._p->_indices->_ids->size()>1)) && !pair.second._p->is_matrix()) {
//                    auto dim = pair.second._p->get_dim(i);
//                    if (pair.second._sign) {
//                        for (size_t j = 0; j<dim; j++) {
//                            res += eval_coef(pair.second._coef,i,j) * ((param<double>*)(pair.second._p.get()))->eval(i,j);
//                        }
//                    }
//                    else {
//                        for (size_t j = 0; j<dim; j++) {
//                            res -= eval_coef(pair.second._coef,i,j) * ((param<double>*)(pair.second._p.get()))->eval(i,j);
//                        }
//                    }
//                }
//                else {
//                    if (pair.second._sign) {
//                        res += eval_coef(pair.second._coef,i) * ((param<double>*)(pair.second._p.get()))->eval(i);
//                    }
//                    else {
//                        res -= eval_coef(pair.second._coef,i) * ((param<double>*)(pair.second._p.get()))->eval(i);
//                    }
//                }
//            }
////        }
//        //                res += eval_lterms(i);
////        if(!_qterms->empty()){
//            for (auto &pair:*_qterms) {
//                double qval = 0;
//                if (pair.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
//                    assert(pair.second._p->first->_dim[1]==1 && pair.second._coef->_dim[0]==pair.second._p->second->_dim[0]);
//                    for (auto i = 0; i<pair.second._p->first->_dim[0]; i++) {
//                        for (auto j = 0; j<pair.second._p->first->_dim[0]; j++) {
//                            qval += eval_coef(pair.second._coef,i,j) * ((param<double>*)(pair.second._p->first.get()))->eval(i) * ((param<double>*)(pair.second._p->second.get()))->eval(j);
//                        }
//                    }
//                }
//                else if (pair.second._p->first->is_matrix() && !pair.second._p->second->is_matrix() && !pair.second._p->second->_is_transposed) {//matrix * vect
//                    for (size_t j = 0; j<pair.second._p->second->_dim[0]; j++) {
//                        qval += ((param<double>*)(pair.second._p->first.get()))->eval(i,j) * ((param<double>*)(pair.second._p->second.get()))->eval(j);
//                    }
//                    qval *= eval_coef(pair.second._coef,i);
//                }
//                else if (!pair.second._p->first->is_matrix() && pair.second._p->first->_is_transposed && pair.second._p->second->is_matrix() ) {//transposed vect * matrix
//                    for (size_t j = 0; j<pair.second._p->first->_dim[0]; j++) {
//                        qval += ((param<double>*)(pair.second._p->first.get()))->eval(j) * ((param<double>*)(pair.second._p->second.get()))->eval(j,i);
//                    }
//                    qval *= eval_coef(pair.second._coef,i);
//                }
//                else if (!pair.second._p->first->is_matrix() && pair.second._p->first->_is_transposed && !pair.second._p->second->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
//                    for (size_t j = 0; j<pair.second._p->first->_dim[1]; j++) {
//                        qval += ((param<double>*)(pair.second._p->first.get()))->eval(j) * ((param<double>*)(pair.second._p->second.get()))->eval(j);
//                    }
//                    qval *= eval_coef(pair.second._coef,i);
//                }
//                else if (!pair.second._coef->is_matrix() && pair.second._coef->_is_transposed && !pair.second._p->first->is_matrix()) {//transposed vect * vec, a dot product of two vectors
//                    for (size_t j = 0; j<pair.second._p->first->_dim[0]; j++) {
//                        qval += eval_coef(pair.second._coef,j) * ((param<double>*)(pair.second._p->first.get()))->eval(j) * ((param<double>*)(pair.second._p->second.get()))->eval(j);
//                    }
//                }
//                else {
//                    qval += eval_coef(pair.second._coef,i) * ((param<double>*)(pair.second._p->first.get()))->eval(i) * ((param<double>*)(pair.second._p->second.get()))->eval(i);
//                }
//                if (!pair.second._sign) {
//                    qval *= -1;
//                }
//                res += qval;
//            }
////        }
//        //                res += eval_qterms(i);
////        if(!_pterms->empty()){
//            for (auto &pair:*_pterms) {
//                double pval = 1;
//                for (auto &vpair: *pair.second._l) {
//                    pval *= std::pow(eval(vpair.first, i), vpair.second);
//                }
//                pval *= eval_coef(pair.second._coef,i);
//                if (!pair.second._sign) {
//                    pval *= -1;
//                }
//                res += pval;
//            }
//            //                res += eval_pterms(i);
////        }
//        if(_expr)
//            res += eval_expr(_expr,i);
//        if (_vars->empty() && _params->empty()) {
//            _val->at(0) = res;
//            _evaluated = true;
//        }
//        else {
//            if (_vars->empty() && i==_val->size()-1) {
//                _evaluated = true;
//            }
//            _val->at(i) = res;
//        }
//        return res;
//    }
//
//
//
//    void func_::delete_var(const string& vid){
//        auto vit = _vars->find(vid);
//        if (vit==_vars->end()) {
//            return;
//        }
//        _vars->erase(vit);
//    }
//
//    void func_::delete_param(const string& vid){
//        auto vit = _params->find(vid);
//        if (vit==_params->end()) {
//            return;
//        }
//        _params->erase(vit);
//    }
//
//
//
//
    unsigned func_::nb_occ_var(string name) const{/**< Returns the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(name);
        if (pair_it==_vars->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }

    unsigned func_::nb_occ_param(string name) const{/**< Returns the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(name);
        if (pair_it==_params->end()) {
            return 0;
        }
        else {
            return get<1>(*pair_it).second;
        }
    }

    void func_::incr_occ_var(string str){/**< Increases the number of occurences the variable has in this function. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }

    void func_::incr_occ_param(string str){/**< Increases the number of occurences the parameter has in this function. */
        auto pair_it = _params->find(str);
        if (pair_it==_params->end()) {
            throw invalid_argument("Non-existing variable in function!\n");
        }
        else {
            get<1>(*pair_it).second++;
        }
    }

    void func_::decr_occ_var(string str, int nb){/**< Decreases the number of occurences the variable has in this function by nb. */
        auto pair_it = _vars->find(str);
        if (pair_it==_vars->end()) {
            return;
        }
        else {
            get<1>(*pair_it).second-=nb;
            if (get<1>(*pair_it).second==0) {
                _vars->erase(pair_it);
            }
        }
    }

    void func_::decr_occ_param(string str, int nb){/**< Decreases the number of occurences the parameter has in this function by nb. */
        auto pair_it = _params->find(str);
        if (pair_it==_params->end()) {
            return;
        }
        else {
            get<1>(*pair_it).second -= nb;
            if (get<1>(*pair_it).second==0) {
                _params->erase(pair_it);
            }
        }
    }


    


    bool func_::is_convex(size_t idx) const{
        return (_convexity->at(idx)==convex_ || _convexity->at(idx)==linear_);
    };

    bool func_::is_concave(size_t idx) const{
        return (_convexity->at(idx)==concave_ || _convexity->at(idx)==linear_);
    };


//    bool func_::is_number() const{
//        return (_vars->empty() && _params->empty());
//    }

    

    bool func_::is_linear() const{
        return (_ftype==lin_);
    };

    bool func_::is_quadratic() const{
        return (_ftype==quad_);
    };

    bool func_::is_polynomial() const{
        return (_ftype==pol_);
    };

    bool func_::is_nonlinear() const{
        return (_ftype==nlin_);
    };

//    bool func_::is_complex() const{
//        for(auto &it: *_vars){
//            if (it.second.first->is_complex()) {
//                return true;
//            }
//        }
//        for(auto &it: *_params){
//            if (it.second.first->is_complex()) {
//                return true;
//            }
//        }
//        return false;
//    };
//
//
//    bool func_::is_zero() const{/*<< A function is zero if it is constant and equals zero or if it is a sum of zero valued parameters */
//        if(is_number() && !_is_vector && !_is_transposed && t_eval(this)==0){
//            return true;
//        }
//        if (_ftype==const_ && _cst->is_zero()){
//            for (auto& it:*_params) {
//                if (!it.second.first->is_zero()) {
//                    return false;
//                }
//            }
//            return true;
//        }
//        return false;
//    }
//
//    bool func_::is_transposed() const {
//        return _is_transposed;
//    }
//
    FType func_::get_ftype() const { return _ftype;}
    
//
//
    qterm* func_::get_square(shared_ptr<param_> p){ /**< Returns the quadratic term containing a square of p or nullptr if none exists. **/
        for (auto pair_it = _qterms->begin(); pair_it != _qterms->end(); pair_it++) {
            if (pair_it->second._p->first==p && pair_it->second._p->second==p) {
                return &pair_it->second;
            }
        }
        return nullptr;
    }
//
//    func_ func_::get_outer_app(){ /**< Returns an outer-approximation of the function using the current value of the variables **/
//        func_ res; // res = gradf(x*)*(x-x*) + f(x*)
//        param_* v;
//        for(auto &it: *_vars){
//            v = it.second.first.get();
//            res += (get_stored_derivative(v->_name)->eval())*((*v) - t_eval(v));
//        }
//        res += eval();
//        return res;
//    }
//
//
//    pair<double,double>* func_::get_all_range() const{
//        return _all_range;
//    }
//
//    Sign func_::get_all_sign() const{
//        return _all_sign;
//    }
//
//    Sign func_::get_sign(size_t idx) const{
//        return _sign->at(idx);
//    }
//
//    Sign func_::get_all_sign(const lterm& l) {
//        if (l._coef->is_zero()) {
//            return zero_;
//        }
//        if (l._coef->get_all_sign()==unknown_ || l._p->get_all_sign()==unknown_) {
//            return unknown_;
//        }
//        auto s = l._coef->get_all_sign() * l._p->get_all_sign();
//        if(s == 1 || s == 2) {
//            if (l._sign) {
//                return non_neg_;
//            }
//            else {
//                return non_pos_;
//            }
//        }
//        if(s == 4) {
//            if (l._sign) {
//                return pos_;
//            }
//            else {
//                return neg_;
//            }
//        }
//        if(s == -1 || s == -2) {
//            if (l._sign) {
//                return non_pos_;
//            }
//            else{
//                return non_neg_;
//            }
//        }
//        if(s == -4) {
//            if (l._sign) {
//                return neg_;
//            }
//            else {
//                return pos_;
//            }
//        }
//        return unknown_;
//    }
//
//    Sign func_::get_all_sign(const qterm& l) {
//        if (l._coef->is_zero()) {
//            return zero_;
//        }
//        if (l._coef->get_all_sign()==unknown_ || l._p->first->get_all_sign()==unknown_ || l._p->second->get_all_sign()==unknown_) {
//            return unknown_;
//        }
//        auto s = l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign();
//        if(s == 1 || s == 2 || s == 4) {
//            if (l._sign) {
//                return non_neg_;
//            }
//            else {
//                return non_pos_;
//            }
//        }
//        if(s == 8) {
//            if (l._sign) {
//                return pos_;
//            }
//            else {
//                return neg_;
//            }
//        }
//        if(s == -1 || s == -2 || s == -4) {
//            if (l._sign) {
//                return non_pos_;
//            }
//            else{
//                return non_neg_;
//            }
//        }
//        if(s == -8) {
//            if (l._sign) {
//                return neg_;
//            }
//            else {
//                return pos_;
//            }
//        }
//        return unknown_;
//    }
//
//    Sign func_::get_all_sign(const pterm& l) {
//        if (l._coef->is_zero()) {
//            return zero_;
//        }
//        //        if (l._coef->get_all_sign()==unknown_ || l._p->first->get_all_sign()==unknown_ || l._p->second->get_all_sign()==unknown_) {
//        //            return unknown_;
//        //        }
//        //        if (l._sign) {
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 2) {
//        //                return non_neg_;
//        //            }
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 4) {
//        //                return pos_;
//        //            }
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -2) {
//        //                return non_pos_;
//        //            }
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -4) {
//        //                return neg_;
//        //            }
//        //        }
//        //        else {
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 2) {
//        //                return non_pos_;
//        //            }
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == 4) {
//        //                return neg_;
//        //            }
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -2) {
//        //                return non_neg_;
//        //            }
//        //            if(l._coef->get_all_sign() * l._p->first->get_all_sign() * l._p->second->get_all_sign() == -4) {
//        //                return pos_;
//        //            }
//        //        }
//        return unknown_;
//    }
//
//
//
    
    
    
    /**
     Index the function and its variables/parameters using nodes of a graph
     @param[in] vec vector of nodes
     @return current function
     */
//    func_& func_::in(const vector<Node*>& vec){
//        _nb_vars = 0;
//        string key;
//        auto iter = _vars->begin();
//        while (iter!=_vars->end()) {
//            auto pair = (*iter++);
//            auto v = pair.second.first;
//            if(!v->_indices){
//                v->index_in(vec);
//            }
//            else if(v->_indices->_type==in_arcs_){
//                v->index_in_arcs(vec);
//            }
//            else if(v->_indices->_type==out_arcs_){
//                v->index_out_arcs(vec);
//            }
//            else if(v->_indices->_type==in_gens_){
//                v->index_in_aux(vec,"gens");
//            }
//            if (!v->_is_vector) {// i.e., it is not transposed
//                _nb_vars++;
//            }
//            else {
//                _nb_vars += v->get_dim();
//            }
//        }
//        iter = _params->begin();
//        while (iter!=_params->end()) {
//            auto pair = (*iter++);
//            auto p = pair.second.first;
//            if(!p->_indices){
//                p->index_in(vec);
//            }
//            else if(p->_indices->_type==in_arcs_){
//                p->index_in_arcs(vec);
//            }
//            else if(p->_indices->_type==out_arcs_){
//                p->index_out_arcs(vec);
//            }
//            else if(p->_indices->_type==in_gens_){
//                p->index_in_aux(vec,"gens");
//            }
//        }
//        return *this;
//    }
    
    /**
     Index the function and its variables/parameters using the indices in ids
     @param[in] ids indices
     @return current function
     */
//    func_& func_::in(const indices& ids) {
//        _nb_vars = 0;
//        string key;
//        auto iter = _vars->begin();
//        while (iter!=_vars->end()) {
//            auto pair = (*iter++);
//            auto v = pair.second.first;
//            if(!v->is_indexed()){
//                v->index_in(ids);
//            }
//            if (!v->_is_vector) {// i.e., it is not transposed
//                _nb_vars++;
//            }
//            else {
//                _nb_vars += v->get_dim();
//            }
//        }
//        iter = _params->begin();
//        while (iter!=_params->end()) {
//            auto pair = (*iter++);
//            auto p = pair.second.first;
//            if(!p->is_indexed()){
//                p->index_in(ids);
//            }
//        }
//        _indices = make_shared<indices>(ids);
//        return *this;
//    }
    
    

        
//
//
//    string poly_to_str(const constant_* c, size_t inst){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//
//        if (!c) {
//            return "null";
//        }
//        switch (c->get_type()) {
//            case binary_c: {
//                return ((constant<bool>*)(c))->to_str();
//                break;
//            }
//            case short_c: {
//                return ((constant<short>*)(c))->to_str();
//                break;
//            }
//            case integer_c: {
//                return ((constant<int>*)(c))->to_str();
//                break;
//            }
//            case float_c: {
//                return ((constant<float>*)(c))->to_str();
//                break;
//            }
//            case double_c: {
//                return ((constant<double>*)(c))->to_str();
//                break;
//            }
//            case long_c: {
//                return ((constant<long double>*)(c))->to_str();
//                break;
//            }
//            case complex_c: {
//                return ((constant<Cpx>*)(c))->to_str();
//                break;
//            }
//            case par_c:{
//                auto p_c = (param_*)(c);
//                switch (p_c->get_intype()) {
//                    case binary_:
//                        return ((param<bool>*)p_c)->to_str(inst);
//                        break;
//                    case short_:
//                        return ((param<short>*)p_c)->to_str(inst);
//                        break;
//                    case integer_:
//                        return ((param<int>*)p_c)->to_str(inst);
//                        break;
//                    case float_:
//                        return ((param<float>*)p_c)->to_str(inst);
//                        break;
//                    case double_:
//                        return ((param<double>*)p_c)->to_str(inst);
//                        break;
//                    case long_:
//                        return ((param<long double>*)p_c)->to_str(inst);
//                        break;
//                    case complex_:
//                        return ((param<Cpx>*)p_c)->to_str(inst);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case uexp_c: {
//                return ((uexpr*)c)->to_str(inst);
//                break;
//            }
//            case bexp_c: {
//                return ((bexpr*)c)->to_str(inst);
//                break;
//            }
//            case var_c: {
//                auto p_c = (param_*)(c);
//                switch (p_c->get_intype()) {
//                    case binary_:
//                        return ((var<bool>*)p_c)->get_name(inst);
//                        break;
//                    case short_:
//                        return ((var<short>*)p_c)->get_name(inst);
//                        break;
//                    case integer_:
//                        return ((var<int>*)p_c)->get_name(inst);
//                        break;
//                    case float_:
//                        return ((var<float>*)p_c)->get_name(inst);
//                        break;
//                    case double_:
//                        return ((var<double>*)p_c)->get_name(inst);
//                        break;
//                    case long_:
//                        return ((var<long double>*)p_c)->get_name(inst);
//                        break;
//                    case complex_:
//                        return ((var<Cpx>*)p_c)->get_name(inst);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case func_c: {
//                return ((func_*)c)->to_str(inst);
//                break;
//            }
//            default:
//                break;
//        }
//        return "null";
//    }
//
//    string poly_to_str(const constant_* c, size_t inst1, size_t inst2){/**< printing c, detecting the right class, i.e., constant<>, param<>, uexpr or bexpr. */
//
//        if (!c) {
//            return "null";
//        }
//        switch (c->get_type()) {
//            case binary_c: {
//                return ((constant<bool>*)(c))->to_str();
//                break;
//            }
//            case short_c: {
//                return ((constant<short>*)(c))->to_str();
//                break;
//            }
//            case integer_c: {
//                return ((constant<int>*)(c))->to_str();
//                break;
//            }
//            case float_c: {
//                return ((constant<float>*)(c))->to_str();
//                break;
//            }
//            case double_c: {
//                return ((constant<double>*)(c))->to_str();
//                break;
//            }
//            case long_c: {
//                return ((constant<long double>*)(c))->to_str();
//                break;
//            }
//            case complex_c: {
//                return ((constant<Cpx>*)(c))->to_str();
//                break;
//            }
//            case par_c:{
//                auto p_c = (param_*)(c);
//                switch (p_c->get_intype()) {
//                    case binary_:
//                        return ((param<bool>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    case short_:
//                        return ((param<short>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    case integer_:
//                        return ((param<int>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    case float_:
//                        return ((param<float>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    case double_:
//                        return ((param<double>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    case long_:
//                        return ((param<long double>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    case complex_:
//                        return ((param<Cpx>*)p_c)->to_str(inst1,inst2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//            case uexp_c: {
//                return ((uexpr*)c)->to_str(inst2);
//                break;
//            }
//            case bexp_c: {
//                return ((bexpr*)c)->to_str(inst2);
//                break;
//            }
//            case var_c: {
//                auto p_c = (param_*)(c);
//                switch (p_c->get_intype()) {
//                    case binary_:
//                        return ((var<bool>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    case short_:
//                        return ((var<short>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    case integer_:
//                        return ((var<int>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    case float_:
//                        return ((var<float>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    case double_:
//                        return ((var<double>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    case long_:
//                        return ((var<long double>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    case complex_:
//                        return ((var<Cpx>*)p_c)->get_name(inst1,inst2);
//                        break;
//                    default:
//                        break;
//                }
//                break;
//            }
//
//            case func_c: {
//                if(c->is_matrix()){
//                    return to_string(((func_*)c)->get_val(inst1,inst2));
//                }
//                return ((func_*)c)->to_str(inst2);
//                break;
//            }
//            default:
//                break;
//        }
//        return "null";
//    }
//
//
//    func_ conj(const func_& f){
//        func_ newf(f);
//        newf._is_conjugate = !newf._is_conjugate;
//        return newf;
//    }
//
//    func_ ang(const func_& f){
//        func_ newf(f);
//        newf._is_angle = true;
//        return newf;
//    }
//
//    func_ sqrmag(const func_& f){
//        func_ newf(f);
//        newf._is_sqrmag = true;
//        return newf;
//    }
//
//    func_ real(const func_& f){
//        func_ newf(f);
//        newf._is_real = true;
//        return newf;
//    }
//
//    func_ imag(const func_& f){
//        func_ newf(f);
//        newf._is_imag = true;
//        return newf;
//    }
//
    
//
//    template func<double> constant<double>::get_real() const;
//    template func<double> constant<double>::get_imag() const;
//    template func<double> constant<Cpx>::get_real() const;
//    template func<double> constant<Cpx>::get_imag() const;
//    template func<double> var<double>::get_real() const;
//    template func<double> var<double>::get_imag() const;
//    template func<double> var<Cpx>::get_real() const;
//    template func<double> var<Cpx>::get_imag() const;
//    template <>
//    func<double> var<double>::get_real() const;
//    template <>
//    func<double> var<double>::get_imag() const;
//    template <>
//    func<double> param<double>::get_real() const;
//    template <>
//    func<double> param<double>::get_imag() const;
//    
//    template <>
//    func<double> constant<Cpx>::get_real() const;
//    template <>
//    func<double> constant<Cpx>::get_imag() const;
//    template <>
//    func<double> var<Cpx>::get_real() const;
//    template <>
//    func<double> var<Cpx>::get_imag() const;
//    template <>
//    func<double> param<Cpx>::get_real() const;
//    template <>
//    func<double> param<Cpx>::get_imag() const;
}
