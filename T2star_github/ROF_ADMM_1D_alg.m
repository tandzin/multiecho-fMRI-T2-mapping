function [u_TV, out]=ROF_ADMM_1D_alg(b,u_0,mu,beta,gamma,delta,nr_shrinks,nr_grads,nr_threads)

mu2=mu/2;
beta2=beta/2;
beta_inv=1/beta;
image_size=size(b);
out.alpha_k = [];
out.TV_term = [];
out.TV_norm = [];
out.w_norm = [];
w_norm=0;
out.data_term = [];
out.Q_kcu = [];
out.gradient_count = [];
out.shrink_count = [];
out.P_k = [];
out.Lambda_k = [];

u_k=u_0;%(1)initialize image
e_k=u_k-b;%(1) initialize image error
[zV_n]=Dv_U(u_0,nr_threads);%(2)D times initial u_0, components of image gradient wrt. x,y

wV_n=zeros(image_size);%(3)components of TV estimate
cV_n=zeros(image_size);%(3) components of scaled Lagrange coefficients for TV mismatch
dV_k=zV_n;%(4)y component of quadratic TV term d_k=z_k-w_k-c_k;

Lambda_k=beta2*(sum(zV_n(:).*zV_n(:)));%initializer

for shrink_count=1:nr_shrinks%denoted with n
    disp(['shrink_count=' num2str(shrink_count)])
    for gradient_count=1:nr_grads
        g_k=beta*DT_Gv(dV_k,nr_threads)+mu*e_k;%(5),sum(d'D'), gradient of Lambda w.r.t. u
        %-----------------------CALCULATE STEP SIZE alpha--------------------------
        %calculate stepsize (7)
        %calculate the denominator of (7)
        [dvg_k]=Dv_U(g_k,nr_threads);%(6)D.g_k in denominator of (7)
        g_k_tv=beta*DT_Gv(dvg_k,nr_threads);%(6)g_k'.D' in denominator of (7)
        clear dvg_k
        g_k_data=mu*g_k;
        R_g_k=g_k_tv+g_k_data;%
        g_k_t_R_g_k=sum(g_k(:).*R_g_k(:));
        clear g_k_tv g_k_data R_g_k
        g_k_square=sum(g_k(:).*g_k(:));%numerator of stepsize(7)
        alpha_k=g_k_square/g_k_t_R_g_k;
        %         disp(['alpha_k = ' num2str(alpha_k)]);
        alpha_k=gamma*alpha_k;
        out.alpha_k = [out.alpha_k; alpha_k];
        %-----------------end of CALCULATE STEP SIZE alpha--------------------------
        u_k=u_k-alpha_k*g_k;%(8)
        [zV_n]=Dv_U(u_k,nr_threads);%(9)D.u_k, x,y-gradient of the new u_k
        TV_norm=sum(sum(abs(zV_n(:))));

        rV_k=zV_n-wV_n;%(10) horizontal TV error 
        dV_k=rV_k-cV_n;%(10)vertical scaled Lagrangian error d_i for all i
        TV_term=beta2*(sum(dV_k(:).*dV_k(:)));%(11)
        e_k=u_k-b;
        e_k_square=sum(e_k(:).*e_k(:));%(12)
        data_term=mu2*e_k_square;
        Q_kcu=TV_term+data_term;%Q value after c and u update
        P_k=TV_norm+data_term;%(16)
        out.TV_norm = [out.TV_norm; TV_norm];
        out.w_norm = [out.w_norm;w_norm];
        out.TV_term = [out.TV_term; TV_term];
        out.data_term = [out.data_term; data_term];
        out.Q_kcu = [out.Q_kcu;Q_kcu];
        out.gradient_count = [out.gradient_count; gradient_count];
        out.shrink_count = [out.shrink_count; shrink_count];
        out.P_k = [out.P_k;P_k];
        out.Lambda_k = [out.Lambda_k;Lambda_k];
    end % for gradient_count, denoted with k, 
    % shrinkage
    sV_n=zV_n-cV_n;%(13)

    wV_n = max(abs(sV_n)-beta_inv,0).*sign(sV_n);%(14)
    w_norm=sum(sum(abs(wV_n(:))));
    clear sV_n
    rV_k=zV_n-wV_n;%(10) residual
    dV_k=rV_k-cV_n;%(10)vertical scaled Lagrangian error d_i for all i
    TV_term=beta2*(sum(dV_k(:).*dV_k(:)));
    Q_kw=TV_term+data_term;%(15), %Q value after w update
    L_k=w_norm+Q_kw;%(17)
    Lambda_k=L_k-beta2*(sum(cV_n(:).*cV_n(:)));%(17)
    cV_n=cV_n-delta*(zV_n-wV_n);%(19) delta slows down the integrator

    dV_k=rV_k-cV_n;%(20)
    TV_term_d=beta2*(sum(dV_k(:).*dV_k(:)));
    Q_k_d=TV_term_d+data_term;
    L_k_d=w_norm+Q_k_d;
    Lambda_k_d=L_k_d-beta2*(sum(cV_n(:).*cV_n(:)));
end %shrink_count, denoted with n
u_TV=u_k;
end % ROF_ADMM_alg

