//function graphic web teplate for runge kutta,but...

void finite_diff_method(const mu::string_type& exp,double x,int i){
    
    mu::Parser parser;
    parser.SetExpr(exp);  
    x +=h;  
    parser.DefineVar(mu::string_type(L"x"),&x);

    points[i].dy_dx = (parser.Eval() - points[i].y) / h;
    
    if(i>0)
    {
        points[i].d2y_dx2 = (points[i].dy_dx - points[i-1].dy_dx) / h;
    }

}