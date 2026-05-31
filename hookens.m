clc; clear all; close all

%% ═══════════════════════════════════════════════════════════
%  PARÁMETROS DEL MECANISMO
%% ═══════════════════════════════════════════════════════════
r1=0.1; r2=0.05; r3=0.25; r4=0.125;

phi20   = 2.5;
omega20 = degtorad(10);
alpha20 = 0;

M1=1; M2=1; M3=1; M4=1;

I1=1/12*M1*r1^2;
I2=1/12*M2*r2^2;
I3=1/12*M3*r3^2;
I4=1/12*M4*r4^2;
g=0;

Mmat = diag([M1 M1 I1 M2 M2 I2 M3 M3 I3 M4 M4 I4]);
Fext = [0;-M1*g;0;  0;-M2*g;0;  0;-M3*g;0;  0;-M4*g;0];

q0 = [1; 0; 0; -0.012; 0.022; 2.0944; 0.025; 0.080; 0.67963121; 0.087; 0.060; 1.79472207];

k  = 2;       % N/m
L0 = 0.1376;  % m

%% ═══════════════════════════════════════════════════════════
%  FASE 1 — FOTO EN t=0
%% ═══════════════════════════════════════════════════════════
q = q0; t = 0;
tol = 100; iter = 0;
while tol > 1e-7 && iter < 100
    iter = iter+1;
    [Fq, J] = ecuaciones(q, t, r1, r2, r3, r4, phi20, omega20, alpha20);
    q   = -inv(J)*Fq + q;
    tol = norm(Fq,2);
    q_foto = q;
end

O2_0 = [q(1);q(2)] + Rot(q(3))*[-r1/2; 0];
O4_0 = [q(1);q(2)] + Rot(q(3))*[ r1/2; 0];
A_0  = [q(4);q(5)] + Rot(q(6))*[ r2/2; 0];
B_0  = [q(10);q(11)] + Rot(q(12))*[r4/2; 0];

d1   = B_0 - A_0;
d2   = O4_0 - O2_0;
Amat = [d1, -d2];
bvec = O2_0 - A_0;
if abs(det(Amat)) < 1e-10
    error('CIV 2-4 no está definido en t=0 (barras paralelas).')
end
params     = Amat \ bvec;
CIV24_foto = A_0 + params(1)*d1;

fprintf('CIV 2-4 en t=0: [%.5f,  %.5f]\n', CIV24_foto(1), CIV24_foto(2))

G2_0 = [q(4);  q(5)];
G4_0 = [q(10); q(11)];

Ps2_0 = CIV24_foto;
Ps4_0 = CIV24_foto;

poi_r2_local = Rot(q(6))'  * (Ps2_0 - G2_0);
poi_r4_local = Rot(q(12))' * (Ps4_0 - G4_0);

fprintf('Distancia inicial entre anclajes: %.5f m  (L0 = %.5f m)\n', ...
        norm(Ps4_0 - Ps2_0), L0)

%% ═══════════════════════════════════════════════════════════
%  FASE 2 — SOLVER COMPLETO
%% ═══════════════════════════════════════════════════════════
q = q0; contador = 0;
 
for t = 0:0.1:36
%for t=0:0.1:60
%for t = 0:0.1:45
    contador = contador+1;
 
    tol = 100; iter = 0;
    while tol > 1e-7 && iter < 100
        iter = iter+1;
        [Fq, J] = ecuaciones(q, t, r1, r2, r3, r4, phi20, omega20, alpha20);
        q   = -inv(J)*Fq + q;
        tol = norm(Fq,2);
    end
    if iter >= 99, fprintf('No convergió en t=%.2f\n', t); end
 
    P(:,contador)        = q;
    time(contador)       = t;
    J_hist(:,:,contador) = J;
 
    v  = -inv(J)*[0;0;0;0;0;0;0;0;0;0;0; -alpha20*t - omega20];
    V(:,contador) = v;
 
    Jp = jacobian_punto(q, v, r1, r2, r3, r4);
    a  = -inv(J)*(Jp*v + [0;0;0;0;0;0;0;0;0;0;0; -alpha20]);
    Ac(:,contador) = a;
 
    G2_i = [q(4);  q(5)];
    G4_i = [q(10); q(11)];
    Ps2  = G2_i + Rot(q(6))  * poi_r2_local;
    Ps4  = G4_i + Rot(q(12)) * poi_r4_local;
 
    L_actual = norm(Ps4 - Ps2);
    delta    = L_actual - L0;
    if L_actual > 1e-10
        u_spring = (Ps4 - Ps2) / L_actual;
    else
        u_spring = [0; 0];
    end
    F_spring = k * delta * u_spring;
 
    dRot2 = [-sin(q(6))  -cos(q(6));  cos(q(6))  -sin(q(6))];
    dRot4 = [-sin(q(12)) -cos(q(12)); cos(q(12)) -sin(q(12))];
 
    Q_spring = zeros(12,1);
    Q_spring(4)  =  F_spring(1);
    Q_spring(5)  =  F_spring(2);
    Q_spring(6)  =  (dRot2 * poi_r2_local)' * F_spring;
    Q_spring(10) = -F_spring(1);
    Q_spring(11) = -F_spring(2);
    Q_spring(12) = -(dRot4 * poi_r4_local)' * F_spring;
 
    lambda           = inv(J') * (Mmat*a - Fext - Q_spring);
    LBA(:,contador)  = lambda;
 
    lambda_base           = inv(J') * (Mmat*a - Fext);
    Torque_base(contador) = lambda_base(12);
 
    elongacion(contador) = delta;
    energia(contador)    = 0.5 * k * delta^2;
    longitud(contador)   = L_actual;
    traj_s2(:,contador)  = Ps2;
    traj_s4(:,contador)  = Ps4;
 
    % ── Extremo libre de la barra 3 ───────────────────────────────────
    traj_tip3(:, contador) = [q(7); q(8)] + Rot(q(9)) * [r3/2; 0];
 
    % CIV 2-4
    O2_i = [q(1);q(2)] + Rot(q(3))*[-r1/2; 0];
    O4_i = [q(1);q(2)] + Rot(q(3))*[ r1/2; 0];
    A_i  = [q(4);q(5)] + Rot(q(6))*[ r2/2; 0];
    B_i  = [q(10);q(11)] + Rot(q(12))*[r4/2; 0];
 
    d1 = B_i - A_i;  d2 = O4_i - O2_i;
    Am = [d1, -d2];  bv = O2_i - A_i;
    if abs(det(Am)) > 1e-10
        par = Am \ bv;
        CIV24(:,contador) = A_i + par(1)*d1;
    elseif contador > 1
        CIV24(:,contador) = CIV24(:,contador-1);
    else
        CIV24(:,contador) = [NaN; NaN];
    end
 
    d3  = A_i - O2_i;  d4  = B_i - O4_i;
    Am2 = [d3, -d4];   bv2 = O4_i - O2_i;
    if abs(det(Am2)) > 1e-10
        par2 = Am2 \ bv2;
        CIV13(:,contador) = O2_i + par2(1)*d3;
    elseif contador > 1
        CIV13(:,contador) = CIV13(:,contador-1);
    else
        CIV13(:,contador) = [NaN; NaN];
    end
end

% % % ═══════════════════════════════════════════════════════════
% % %  FASE 3 — VARIABLES DERIVADAS PARA ANÁLISIS
% % % ═══════════════════════════════════════════════════════════
phi2_vec = P(6,:);
N_pasos  = size(P,2);
Ps2_all  = traj_s2;
Ps4_all  = traj_s4;
L_vec    = longitud;
L_min    = min(L_vec);
L_max    = max(L_vec);

p_fixed  = CIV24(:,1);
G2_foto  = [q_foto(4);  q_foto(5)];
G4_foto  = [q_foto(10); q_foto(11)];
r2_loc   = Rot(q_foto(6))'  * (p_fixed - G2_foto);
r4_loc   = Rot(q_foto(12))' * (p_fixed - G4_foto);

N_L0     = 20;
L0_sweep = linspace(0, L_max, N_L0);
colores  = parula(N_L0);

dphi2_v  = diff(phi2_vec);
mask_v   = abs(dphi2_v) > 1e-12;

%── Calcular W vs L0 para el punto CIV24(t=0) ────────────────────────
W_vec     = zeros(1, N_L0);
W_tension = zeros(1, N_L0);
W_compres = zeros(1, N_L0);
TT_store  = zeros(N_L0, N_pasos);   % torque total guardado para graficar

for iL = 1:N_L0
    L0_i         = L0_sweep(iL);
    delta_torque = zeros(1, N_pasos);
    for i = 1:N_pasos
        Li = L_vec(i);
        if Li > 1e-10, u_sp = (Ps4_all(:,i)-Ps2_all(:,i))/Li; else, u_sp=[0;0]; end
        F_sp = k*(Li-L0_i)*u_sp;
        dR2 = [-sin(P(6,i)) -cos(P(6,i)); cos(P(6,i)) -sin(P(6,i))];
        dR4 = [-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
        Q_new=zeros(12,1);
        Q_new(4)=F_sp(1); Q_new(5)=F_sp(2);
        Q_new(6)=(dR2*r2_loc)'*F_sp;
        Q_new(10)=-F_sp(1); Q_new(11)=-F_sp(2);
        Q_new(12)=-(dR4*r4_loc)'*F_sp;
        contrib=J_hist(:,:,i)'\Q_new;
        delta_torque(i)=-contrib(12);
    end
    Torque_total = Torque_base + delta_torque;
    TT_store(iL,:) = Torque_total;

    dT    = diff(Torque_total);
    K_eff = zeros(1,N_pasos-1);
    K_eff(mask_v) = dT(mask_v)./dphi2_v(mask_v);
    neg_m = K_eff < 0;

    delta_vec        = L_vec - L0_i;
    mask_Kneg_t      = neg_m & (delta_vec(1:end-1) > 0);
    mask_Kneg_c      = neg_m & (delta_vec(1:end-1) < 0);

    W_vec(iL)     = sum(abs(dphi2_v(neg_m)));
    W_tension(iL) = sum(abs(dphi2_v(mask_Kneg_t)));
    W_compres(iL) = sum(abs(dphi2_v(mask_Kneg_c)));
end
% % % 
% % % ── Separación horizontal y vertical ─────────────────────────────────
N_sep       = 10;
sep_vec     = (0:N_sep-1)*0.010;
colores_sep = parula(N_sep);

W_sep   = zeros(1,N_sep);
W_sep_v = zeros(1,N_sep);
TT_sep  = zeros(N_sep, N_pasos);
TT_sepv = zeros(N_sep, N_pasos);

for is = 1:N_sep
    d = sep_vec(is);
    %Horizontal
    r2s = Rot(q_foto(6))'  * (CIV24_foto+[d/2;0] - G2_foto);
    r4s = Rot(q_foto(12))' * (CIV24_foto-[d/2;0] - G4_foto);
    %Vertical
    r2v = Rot(q_foto(6))'  * (CIV24_foto+[0;d/2] - G2_foto);
    r4v = Rot(q_foto(12))' * (CIV24_foto-[0;d/2] - G4_foto);

    dts = zeros(1,N_pasos);  dtv = zeros(1,N_pasos);
    for i = 1:N_pasos
        ps2s=[P(4,i);P(5,i)]+Rot(P(6,i))*r2s;  ps4s=[P(10,i);P(11,i)]+Rot(P(12,i))*r4s;
        ps2v=[P(4,i);P(5,i)]+Rot(P(6,i))*r2v;  ps4v=[P(10,i);P(11,i)]+Rot(P(12,i))*r4v;
        Ls=norm(ps4s-ps2s); Lv=norm(ps4v-ps2v);
        if Ls>1e-10, us=(ps4s-ps2s)/Ls; else, us=[0;0]; end
        if Lv>1e-10, uv=(ps4v-ps2v)/Lv; else, uv=[0;0]; end
        Fs=k*(Ls-L0)*us; Fv=k*(Lv-L0)*uv;
        dR2=[-sin(P(6,i)) -cos(P(6,i)); cos(P(6,i)) -sin(P(6,i))];
        dR4=[-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
        Qs=zeros(12,1); Qv=zeros(12,1);
        Qs(4)=Fs(1);Qs(5)=Fs(2);Qs(6)=(dR2*r2s)'*Fs;Qs(10)=-Fs(1);Qs(11)=-Fs(2);Qs(12)=-(dR4*r4s)'*Fs;
        Qv(4)=Fv(1);Qv(5)=Fv(2);Qv(6)=(dR2*r2v)'*Fv;Qv(10)=-Fv(1);Qv(11)=-Fv(2);Qv(12)=-(dR4*r4v)'*Fv;
        cs=J_hist(:,:,i)'\Qs; cv=J_hist(:,:,i)'\Qv;
        dts(i)=-cs(12); dtv(i)=-cv(12);
    end
    TT_sep(is,:)  = Torque_base + dts;
    TT_sepv(is,:) = Torque_base + dtv;

    Ks=zeros(1,N_pasos-1); Kv=zeros(1,N_pasos-1);
    tmp_s=diff(TT_sep(is,:));  Ks(mask_v)=tmp_s(mask_v)./dphi2_v(mask_v);
    tmp_v=diff(TT_sepv(is,:)); Kv(mask_v)=tmp_v(mask_v)./dphi2_v(mask_v);
    W_sep(is)  =sum(abs(dphi2_v(Ks<0)));
    W_sep_v(is)=sum(abs(dphi2_v(Kv<0)));
end
% % % 
% % % ── Mapa de calor ─────────────────────────────────────────────────────
N_grid   = 1504;
N_L0_map = 60;

x_grid = linspace(-0.6, 0.2, N_grid);
y_grid = linspace(-0.3, 0.3, N_grid);

mapa_W_tension = zeros(N_grid,N_grid); mapa_W_compres = zeros(N_grid,N_grid);
mapa_L0_t      = zeros(N_grid,N_grid); mapa_L0_c      = zeros(N_grid,N_grid);
mapa_Lmax      = zeros(N_grid,N_grid); mapa_efic_t    = zeros(N_grid,N_grid);
mapa_efic_c    = zeros(N_grid,N_grid); mapa_L0pct_t   = zeros(N_grid,N_grid);
mapa_L0pct_c   = zeros(N_grid,N_grid);

fprintf('Calculando mapa %dx%d...\n', N_grid, N_grid)

for ix = 1:N_grid
    for iy = 1:N_grid
        p        = [x_grid(ix); y_grid(iy)];
        r2_loc_g = Rot(q_foto(6))'  * (p - G2_foto);
        r4_loc_g = Rot(q_foto(12))' * (p - G4_foto);

        L_g = zeros(1,N_pasos); Ps2_g=zeros(2,N_pasos); Ps4_g=zeros(2,N_pasos);
        for i = 1:N_pasos
            ps2=([P(4,i);P(5,i)]+Rot(P(6,i))*r2_loc_g);
            ps4=([P(10,i);P(11,i)]+Rot(P(12,i))*r4_loc_g);
            Ps2_g(:,i)=ps2; Ps4_g(:,i)=ps4; L_g(i)=norm(ps4-ps2);
        end
        L_max_g=max(L_g);
        if L_max_g<1e-10||L_max_g>2.0, continue; end
        mapa_Lmax(iy,ix)=L_max_g;

        W_best_t=0; L0_best_t=0; W_best_c=0; L0_best_c=0;
        for iL = 1:N_L0_map
            L0_g = L_max_g*(iL-1)/(N_L0_map-1);
            dt_map=zeros(1,N_pasos);
            for i=1:N_pasos
                Li=L_g(i);
                if Li>1e-10, u_sp=(Ps4_g(:,i)-Ps2_g(:,i))/Li; else, u_sp=[0;0]; end
                F_sp=k*(Li-L0_g)*u_sp;
                dR2=[-sin(P(6,i)) -cos(P(6,i)); cos(P(6,i)) -sin(P(6,i))];
                dR4=[-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
                Q_new=zeros(12,1);
                Q_new(4)=F_sp(1);Q_new(5)=F_sp(2);Q_new(6)=(dR2*r2_loc_g)'*F_sp;
                Q_new(10)=-F_sp(1);Q_new(11)=-F_sp(2);Q_new(12)=-(dR4*r4_loc_g)'*F_sp;
                c=J_hist(:,:,i)'\Q_new; dt_map(i)=-c(12);
            end
            TT=Torque_base+dt_map;
            dT=diff(TT); Ke=zeros(1,N_pasos-1); Ke(mask_v)=dT(mask_v)./dphi2_v(mask_v);
            neg_m=Ke<0;
            region_t=calcular_region_eq(TT,Ke,neg_m,L_g,L0_g,k,N_pasos,'tension');
            region_c = calcular_region_eq(TT,Ke,neg_m,L_g,L0_g,k,N_pasos,'compresion');
            dv=L_g-L0_g;
            W_t=sum(abs(dphi2_v(region_t&(dv(1:end-1)>=0))));
            W_c=sum(abs(dphi2_v(region_c&(dv(1:end-1)<0))));
            if W_t>W_best_t, W_best_t=W_t; L0_best_t=L0_g; end
            if W_c>W_best_c, W_best_c=W_c; L0_best_c=L0_g; end
        end
        mapa_W_tension(iy,ix)=W_best_t; mapa_W_compres(iy,ix)=W_best_c;
        mapa_L0_t(iy,ix)=L0_best_t;     mapa_L0_c(iy,ix)=L0_best_c;
        if L_max_g>1e-10
            mapa_efic_t(iy,ix)=W_best_t/L_max_g; mapa_efic_c(iy,ix)=W_best_c/L_max_g;
            mapa_L0pct_t(iy,ix)=100*L0_best_t/L_max_g;
            mapa_L0pct_c(iy,ix)=100*L0_best_c/L_max_g;
        end
    end
    if mod(ix,10)==0, fprintf('  columna %d / %d\n',ix,N_grid); end
end
fprintf('Listo.\n')

[W_t_max,   idx_t]  = max(mapa_W_tension(:)); [iy_t,ix_t]   = ind2sub([N_grid N_grid],idx_t);
[W_c_max,   idx_c]  = max(mapa_W_compres(:)); [iy_c,ix_c]   = ind2sub([N_grid N_grid],idx_c);
[efic_t_max,idx_et] = max(mapa_efic_t(:));    [iy_et,ix_et] = ind2sub([N_grid N_grid],idx_et);
[efic_c_max,idx_ec] = max(mapa_efic_c(:));    [iy_ec,ix_ec] = ind2sub([N_grid N_grid],idx_ec);
% % % 
casos(1).nombre='W_{tens} máximo';    casos(1).ix=ix_t;  casos(1).iy=iy_t;  casos(1).L0=mapa_L0_t(iy_t,ix_t);
casos(2).nombre='W_{comp} máximo';    casos(2).ix=ix_c;  casos(2).iy=iy_c;  casos(2).L0=mapa_L0_c(iy_c,ix_c);
casos(3).nombre='Ef. tensión máxima'; casos(3).ix=ix_et; casos(3).iy=iy_et; casos(3).L0=mapa_L0_t(iy_et,ix_et);
casos(4).nombre='Ef. comp. máxima';   casos(4).ix=ix_ec; casos(4).iy=iy_ec; casos(4).L0=mapa_L0_c(iy_ec,ix_ec);

TC = zeros(4,N_pasos); EK=zeros(4,N_pasos); DK=zeros(4,N_pasos);
region_casos = false(4,N_pasos-1);
phi_eq_casos = nan(1,4);
mask_t_casos = false(4,N_pasos-1);
mask_c_casos = false(4,N_pasos-1);

for ic = 1:4
    p_k  = [x_grid(casos(ic).ix); y_grid(casos(ic).iy)];
    L0_k = casos(ic).L0;
    r2_k = Rot(q_foto(6))'  * (p_k - G2_foto);
    r4_k = Rot(q_foto(12))' * (p_k + [0.0005;0] - G4_foto);
    L_k=zeros(1,N_pasos); dtk=zeros(1,N_pasos);
    for i=1:N_pasos
        ps2=[P(4,i);P(5,i)]+Rot(P(6,i))*r2_k;
        ps4=[P(10,i);P(11,i)]+Rot(P(12,i))*r4_k;
        Li=norm(ps4-ps2); L_k(i)=Li;
        if Li>1e-10, u_sp=(ps4-ps2)/Li; else, u_sp=[0;0]; end
        Fsp=k*(Li-L0_k)*u_sp;
        dR2=[-sin(P(6,i)) -cos(P(6,i)); cos(P(6,i)) -sin(P(6,i))];
        dR4=[-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
        Q=zeros(12,1); Q(4)=Fsp(1);Q(5)=Fsp(2);Q(6)=(dR2*r2_k)'*Fsp;
        Q(10)=-Fsp(1);Q(11)=-Fsp(2);Q(12)=-(dR4*r4_k)'*Fsp;
        c=J_hist(:,:,i)'\Q; dtk(i)=-c(12);
    end
    TC(ic,:)=Torque_base+dtk;
    EK(ic,:)=0.5*k*(L_k-L0_k).^2;
    DK(ic,:)=L_k-L0_k;
    Ke=zeros(1,N_pasos-1); tmp_k=diff(TC(ic,:)); Ke(mask_v)=tmp_k(mask_v)./dphi2_v(mask_v);
    neg_ic=Ke<0;
    sign_ch=find(diff(sign(TC(ic,:))));
    if ~isempty(sign_ch)
        E_cr=zeros(1,numel(sign_ch)); Kn_cr=false(1,numel(sign_ch));
        for sc=1:numel(sign_ch)
            idx_sc=sign_ch(sc);
            E_cr(sc)=(EK(ic,idx_sc)+EK(ic,min(idx_sc+1,N_pasos)))/2;
            if idx_sc<=N_pasos-1, Kn_cr(sc)=neg_ic(idx_sc); end
        end
        E_cd=E_cr; if any(Kn_cr), E_cd(~Kn_cr)=-Inf; end
        [~,bc]=max(E_cd); eq_i=sign_ch(bc);
        phi_eq_casos(ic)=phi2_vec(eq_i)+dphi2_v(min(eq_i,N_pasos-1))* ...
            abs(TC(ic,eq_i))/max(abs(TC(ic,eq_i))+abs(TC(ic,eq_i+1)),1e-12);
    end
    mask_t_casos(ic,:) = calcular_region_eq(TC(ic,:),Ke,neg_ic,L_k,L0_k,k,N_pasos,'tension')    & (DK(ic,1:end-1)>=0);
    mask_c_casos(ic,:) = calcular_region_eq(TC(ic,:),Ke,neg_ic,L_k,L0_k,k,N_pasos,'compresion') & (DK(ic,1:end-1)<0);
    region_casos(ic,:) = mask_t_casos(ic,:) | mask_c_casos(ic,:);
end

%% Puntos óptimos de resorte y separación óptima

%% ═══════════════════════════════════════════════════════════
%  ANÁLISIS DE SEPARACIÓN PARA PUNTOS ÓPTIMOS DEL MAPA
%% ═══════════════════════════════════════════════════════════
% Puntos de interés (compresión y eficiencia en compresión)
% Puedes agregar los de tensión si lo deseas.
p_opt_comp    = [x_grid(ix_c);  y_grid(iy_c)];   % W_comp max
L0_opt_comp   = mapa_L0_c(iy_c, ix_c);
p_opt_comp_ef = [x_grid(ix_ec); y_grid(iy_ec)];  % Ef_comp max
L0_opt_comp_ef= mapa_L0_c(iy_ec, ix_ec);
reg_hor_comp = false(N_sep, N_pasos-1);
delta_hor_comp = zeros(N_sep, N_pasos);

% Vector de separaciones
N_sep       = 10;
sep_vec     = linspace(0, 0.05, N_sep);   % de 0 a 5 cm, ajustable
colores_sep = parula(N_sep);

% Prealocación de resultados
W_hor_comp     = zeros(1, N_sep);
W_ver_comp     = zeros(1, N_sep);
W_hor_comp_ef  = zeros(1, N_sep);
W_ver_comp_ef  = zeros(1, N_sep);

% Opcional: guardar torques totales para gráficas detalladas
TT_hor_comp     = zeros(N_sep, N_pasos);
TT_ver_comp     = zeros(N_sep, N_pasos);
TT_hor_comp_ef  = zeros(N_sep, N_pasos);
TT_ver_comp_ef  = zeros(N_sep, N_pasos);

fprintf('\n===== ANÁLISIS DE SEPARACIÓN =====\n');

for is = 1:N_sep
    d = sep_vec(is);
    fprintf('   Separación d = %.4f m (paso %d/%d)\n', d, is, N_sep);
    
    % --- Caso 1: W_comp max ---
    % Anclajes iniciales para separación horizontal
    p2_h = p_opt_comp + [d/2; 0];
    p4_h = p_opt_comp - [d/2; 0];
    % Anclajes iniciales para separación vertical
    p2_v = p_opt_comp + [0; d/2];
    p4_v = p_opt_comp - [0; d/2];
    
    % Coordenadas locales constantes
    r2_h_loc = Rot(q_foto(6))'  * (p2_h - G2_foto);
    r4_h_loc = Rot(q_foto(12))' * (p4_h - G4_foto);
    r2_v_loc = Rot(q_foto(6))'  * (p2_v - G2_foto);
    r4_v_loc = Rot(q_foto(12))' * (p4_v - G4_foto);
    
    % Contribuciones al par motor
    dt_h = zeros(1, N_pasos);
    dt_v = zeros(1, N_pasos);
    L_h  = zeros(1, N_pasos);
    L_v  = zeros(1, N_pasos);
    for i = 1:N_pasos
        % Posiciones mundiales de los anclajes
        ps2_h = [P(4,i); P(5,i)] + Rot(P(6,i))  * r2_h_loc;
        ps4_h = [P(10,i); P(11,i)] + Rot(P(12,i)) * r4_h_loc;
        ps2_v = [P(4,i); P(5,i)] + Rot(P(6,i))  * r2_v_loc;
        ps4_v = [P(10,i); P(11,i)] + Rot(P(12,i)) * r4_v_loc;
        Lh = norm(ps4_h - ps2_h);
        Lv = norm(ps4_v - ps2_v);
        L_h(i) = Lh;  L_v(i) = Lv;
        
        % Vector unitario y fuerza del resorte (k=2 N/m)
        if Lh > 1e-10, uh = (ps4_h - ps2_h)/Lh; else, uh = [0;0]; end
        if Lv > 1e-10, uv = (ps4_v - ps2_v)/Lv; else, uv = [0;0]; end
        Fh = 2 * (Lh - L0_opt_comp) * uh;
        Fv = 2 * (Lv - L0_opt_comp) * uv;
        
        dR2 = [-sin(P(6,i))  -cos(P(6,i));  cos(P(6,i))  -sin(P(6,i))];
        dR4 = [-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
        
        Qh = zeros(12,1);
        Qh(4)=Fh(1); Qh(5)=Fh(2); Qh(6)=(dR2*r2_h_loc)'*Fh;
        Qh(10)=-Fh(1); Qh(11)=-Fh(2); Qh(12)=-(dR4*r4_h_loc)'*Fh;
        Qv = zeros(12,1);
        Qv(4)=Fv(1); Qv(5)=Fv(2); Qv(6)=(dR2*r2_v_loc)'*Fv;
        Qv(10)=-Fv(1); Qv(11)=-Fv(2); Qv(12)=-(dR4*r4_v_loc)'*Fv;
        
        ch = J_hist(:,:,i)' \ Qh;
        cv = J_hist(:,:,i)' \ Qv;
        dt_h(i) = -ch(12);
        dt_v(i) = -cv(12);
    end
    TT_h = Torque_base + dt_h;
    TT_v = Torque_base + dt_v;
    
    % Cálculo de rigidez efectiva
    dT_h = diff(TT_h);  Ke_h = zeros(1,N_pasos-1); Ke_h(mask_v)=dT_h(mask_v)./dphi2_v(mask_v);
    dT_v = diff(TT_v);  Ke_v = zeros(1,N_pasos-1); Ke_v(mask_v)=dT_v(mask_v)./dphi2_v(mask_v);
    
    % Regiones de rigidez negativa en modo compresión
    reg_h = calcular_region_eq(TT_h, Ke_h, Ke_h<0, L_h, L0_opt_comp, 2, N_pasos, 'compresion');
    reg_v = calcular_region_eq(TT_v, Ke_v, Ke_v<0, L_v, L0_opt_comp, 2, N_pasos, 'compresion');
    
    % Ancho W (suma de Δφ2 en la región, donde además δ<0)
    delta_h = L_h - L0_opt_comp;
    delta_v = L_v - L0_opt_comp;
    W_hor_comp(is) = sum(abs(dphi2_v(reg_h & (delta_h(1:end-1)<0))));
    W_ver_comp(is) = sum(abs(dphi2_v(reg_v & (delta_v(1:end-1)<0))));
    
    TT_hor_comp(is,:) = TT_h;
    TT_ver_comp(is,:) = TT_v;
    
    % --- Caso 2: Ef_comp max (eficiencia) ---
    p2_h_ef = p_opt_comp_ef + [d/2; 0];
    p4_h_ef = p_opt_comp_ef - [d/2; 0];
    p2_v_ef = p_opt_comp_ef + [0; d/2];
    p4_v_ef = p_opt_comp_ef - [0; d/2];
    r2_h_ef_loc = Rot(q_foto(6))'  * (p2_h_ef - G2_foto);
    r4_h_ef_loc = Rot(q_foto(12))' * (p4_h_ef - G4_foto);
    r2_v_ef_loc = Rot(q_foto(6))'  * (p2_v_ef - G2_foto);
    r4_v_ef_loc = Rot(q_foto(12))' * (p4_v_ef - G4_foto);
    
    dt_h_ef = zeros(1, N_pasos);
    dt_v_ef = zeros(1, N_pasos);
    L_h_ef  = zeros(1, N_pasos);
    L_v_ef  = zeros(1, N_pasos);
    for i = 1:N_pasos
        ps2_h = [P(4,i); P(5,i)] + Rot(P(6,i))  * r2_h_ef_loc;
        ps4_h = [P(10,i); P(11,i)] + Rot(P(12,i)) * r4_h_ef_loc;
        ps2_v = [P(4,i); P(5,i)] + Rot(P(6,i))  * r2_v_ef_loc;
        ps4_v = [P(10,i); P(11,i)] + Rot(P(12,i)) * r4_v_ef_loc;
        Lh_ef = norm(ps4_h - ps2_h);
        Lv_ef = norm(ps4_v - ps2_v);
        L_h_ef(i) = Lh_ef;  L_v_ef(i) = Lv_ef;
        
        if Lh_ef > 1e-10, uh = (ps4_h - ps2_h)/Lh_ef; else, uh = [0;0]; end
        if Lv_ef > 1e-10, uv = (ps4_v - ps2_v)/Lv_ef; else, uv = [0;0]; end
        Fh = 2 * (Lh_ef - L0_opt_comp_ef) * uh;
        Fv = 2 * (Lv_ef - L0_opt_comp_ef) * uv;
        
        dR2 = [-sin(P(6,i))  -cos(P(6,i));  cos(P(6,i))  -sin(P(6,i))];
        dR4 = [-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
        
        Qh = zeros(12,1);
        Qh(4)=Fh(1); Qh(5)=Fh(2); Qh(6)=(dR2*r2_h_ef_loc)'*Fh;
        Qh(10)=-Fh(1); Qh(11)=-Fh(2); Qh(12)=-(dR4*r4_h_ef_loc)'*Fh;
        Qv = zeros(12,1);
        Qv(4)=Fv(1); Qv(5)=Fv(2); Qv(6)=(dR2*r2_v_ef_loc)'*Fv;
        Qv(10)=-Fv(1); Qv(11)=-Fv(2); Qv(12)=-(dR4*r4_v_ef_loc)'*Fv;
        
        ch = J_hist(:,:,i)' \ Qh;
        cv = J_hist(:,:,i)' \ Qv;
        dt_h_ef(i) = -ch(12);
        dt_v_ef(i) = -cv(12);
    end
    TT_h_ef = Torque_base + dt_h_ef;
    TT_v_ef = Torque_base + dt_v_ef;
    
    dT_h_ef = diff(TT_h_ef); Ke_h_ef = zeros(1,N_pasos-1); Ke_h_ef(mask_v)=dT_h_ef(mask_v)./dphi2_v(mask_v);
    dT_v_ef = diff(TT_v_ef); Ke_v_ef = zeros(1,N_pasos-1); Ke_v_ef(mask_v)=dT_v_ef(mask_v)./dphi2_v(mask_v);
    
    reg_h_ef = calcular_region_eq(TT_h_ef, Ke_h_ef, Ke_h_ef<0, L_h_ef, L0_opt_comp_ef, 2, N_pasos, 'compresion');
    reg_v_ef = calcular_region_eq(TT_v_ef, Ke_v_ef, Ke_v_ef<0, L_v_ef, L0_opt_comp_ef, 2, N_pasos, 'compresion');
    
    delta_h_ef = L_h_ef - L0_opt_comp_ef;
    delta_v_ef = L_v_ef - L0_opt_comp_ef;
    W_hor_comp_ef(is) = sum(abs(dphi2_v(reg_h_ef & (delta_h_ef(1:end-1)<0))));
    W_ver_comp_ef(is) = sum(abs(dphi2_v(reg_v_ef & (delta_v_ef(1:end-1)<0))));
    
    TT_hor_comp_ef(is,:) = TT_h_ef;
    TT_ver_comp_ef(is,:) = TT_v_ef;


    reg_hor_comp(is,:) = reg_h;
    delta_hor_comp(is,:) = delta_h;
    reg_ver_comp(is,:) = reg_v;
    delta_ver_comp(is,:) = delta_v;

    reg_hor_comp_ef(is,:) = reg_h_ef;
    delta_hor_comp_ef(is,:) = delta_h_ef;
    reg_ver_comp_ef(is,:) = reg_v_ef;
    delta_ver_comp_ef(is,:) = delta_v_ef;
end



% % % 
% % % ── Análisis poloida ──────────────────────────────────────────────────
N_L0_pol=60;
N_pol=N_pasos;
W_best_t_pol=zeros(1,N_pol); W_best_c_pol=zeros(1,N_pol);
Ef_best_t_pol=zeros(1,N_pol); Ef_best_c_pol=zeros(1,N_pol);
L0_opt_t_pol=zeros(1,N_pol); L0_opt_c_pol=zeros(1,N_pol);
Lmax_pol=zeros(1,N_pol); valid_pol=false(1,N_pol);

fprintf('Calculando poloida (%d puntos)...\n', N_pol)
for ip = 1:N_pol
    p_pol=CIV24(:,ip);
    if any(isnan(p_pol))||any(isinf(p_pol)), continue; end
    G2_ip=[P(4,ip);P(5,ip)]; G4_ip=[P(10,ip);P(11,ip)];
    r2lp=Rot(P(6,ip))'*(p_pol-G2_ip); r4lp=Rot(P(12,ip))'*(p_pol-G4_ip);
    L_p=zeros(1,N_pasos); Ps2p=zeros(2,N_pasos); Ps4p=zeros(2,N_pasos);
    for i=1:N_pasos
        ps2=[P(4,i);P(5,i)]+Rot(P(6,i))*r2lp;
        ps4=[P(10,i);P(11,i)]+Rot(P(12,i))*r4lp;
        Ps2p(:,i)=ps2; Ps4p(:,i)=ps4; L_p(i)=norm(ps4-ps2);
    end
    Lmx=max(L_p);
    if Lmx<1e-10||Lmx>2.0, continue; end
    Lmax_pol(ip)=Lmx; valid_pol(ip)=true;
    for iL=1:N_L0_pol
        L0p=Lmx*(iL-1)/(N_L0_pol-1);
        dtp=zeros(1,N_pasos);
        for i=1:N_pasos
            Li=L_p(i);
            if Li>1e-10, u_sp=(Ps4p(:,i)-Ps2p(:,i))/Li; else, u_sp=[0;0]; end
            F_sp=k*(Li-L0p)*u_sp;
            dR2=[-sin(P(6,i)) -cos(P(6,i)); cos(P(6,i)) -sin(P(6,i))];
            dR4=[-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
            Q=zeros(12,1); Q(4)=F_sp(1);Q(5)=F_sp(2);Q(6)=(dR2*r2lp)'*F_sp;
            Q(10)=-F_sp(1);Q(11)=-F_sp(2);Q(12)=-(dR4*r4lp)'*F_sp;
            c=J_hist(:,:,i)'\Q; dtp(i)=-c(12);
        end
        Tp=Torque_base+dtp;
        Ke=zeros(1,N_pasos-1); tmp_p=diff(Tp); Ke(mask_v)=tmp_p(mask_v)./dphi2_v(mask_v);
        neg_p=Ke<0;

        dp      = L_p - L0p;
        reg_t   = calcular_region_eq(Tp,Ke,neg_p,L_p,L0p,k,N_pasos,'tension');
        reg_c   = calcular_region_eq(Tp,Ke,neg_p,L_p,L0p,k,N_pasos,'compresion');
        Wt      = sum(abs(dphi2_v(reg_t & (dp(1:end-1)>=0))));
        Wc      = sum(abs(dphi2_v(reg_c & (dp(1:end-1)<0))));

        if Wt>W_best_t_pol(ip), W_best_t_pol(ip)=Wt; Ef_best_t_pol(ip)=Wt/Lmx; L0_opt_t_pol(ip)=L0p; end
        if Wc>W_best_c_pol(ip), W_best_c_pol(ip)=Wc; Ef_best_c_pol(ip)=Wc/Lmx; L0_opt_c_pol(ip)=L0p; end
    end
    if mod(ip,50)==0, fprintf('  %d / %d\n',ip,N_pol); end
end
fprintf('Listo.\n')

idx_v=find(valid_pol);
x_pol=CIV24(1,idx_v); y_pol=CIV24(2,idx_v);
Wt_v=W_best_t_pol(idx_v); Wc_v=W_best_c_pol(idx_v);
Eft_v=Ef_best_t_pol(idx_v); Efc_v=Ef_best_c_pol(idx_v);
L0t_v=L0_opt_t_pol(idx_v); L0c_v=L0_opt_c_pol(idx_v);
[Wt_max,ip_Wt]=max(Wt_v); [Wc_max,ip_Wc]=max(Wc_v);
[Eft_max,ip_Eft]=max(Eft_v); [Efc_max,ip_Efc]=max(Efc_v);

save("Hookens.mat");
% % % 
% % % ═══════════════════════════════════════════════════════════
% % %  TABLAS EN CONSOLA
% % % ═══════════════════════════════════════════════════════════
fprintf('\n%s\n',repmat('=',1,75))
fprintf('TABLA TENSIÓN — Top 20 puntos del mapa (por eficiencia)\n')
fprintf('%s\n',repmat('=',1,75))
fprintf('%-10s %-10s %-12s %-14s %-14s\n','x [m]','y [m]','L0 [m]','W_tens [rad]','Ef [rad/m]')
fprintf('%s\n',repmat('-',1,75))
[rows_t,cols_t]=find(mapa_W_tension>0);
n_t=length(rows_t); tab_t=zeros(n_t,5);
for ii=1:n_t, r=rows_t(ii);c=cols_t(ii);
    tab_t(ii,:)=[x_grid(c),y_grid(r),mapa_L0_t(r,c),mapa_W_tension(r,c),mapa_efic_t(r,c)]; end
[~,ord_t]=sort(tab_t(:,5),'descend'); tab_t=tab_t(ord_t,:);
for ii=1:min(20,n_t)
    fprintf('%-10.4f %-10.4f %-12.5f %-14.4f %-14.4f\n',tab_t(ii,1),tab_t(ii,2),tab_t(ii,3),tab_t(ii,4),tab_t(ii,5)); end

fprintf('\n%s\n',repmat('=',1,75))
fprintf('TABLA COMPRESIÓN — Top 20 puntos del mapa (por eficiencia)\n')
fprintf('%s\n',repmat('=',1,75))
fprintf('%-10s %-10s %-12s %-14s %-14s\n','x [m]','y [m]','L0 [m]','W_comp [rad]','Ef [rad/m]')
fprintf('%s\n',repmat('-',1,75))
[rows_c,cols_c]=find(mapa_W_compres>0);
n_c=length(rows_c); tab_c=zeros(n_c,5);
for ii=1:n_c, r=rows_c(ii);c=cols_c(ii);
    tab_c(ii,:)=[x_grid(c),y_grid(r),mapa_L0_c(r,c),mapa_W_compres(r,c),mapa_efic_c(r,c)]; end
[~,ord_c]=sort(tab_c(:,5),'descend'); tab_c=tab_c(ord_c,:);
for ii=1:min(20,n_c)
    fprintf('%-10.4f %-10.4f %-12.5f %-14.4f %-14.4f\n',tab_c(ii,1),tab_c(ii,2),tab_c(ii,3),tab_c(ii,4),tab_c(ii,5)); end

fprintf('\n%s\n',repmat('=',1,72))
fprintf('MEJORES PUNTOS DE LA POLOIDA CIV 2-4\n')
fprintf('%s\n',repmat('=',1,72))
fprintf('%-28s  x[m]     y[m]     L0[m]      Valor\n','Criterio')
fprintf('%s\n',repmat('-',1,72))
fprintf('%-28s  %7.4f  %7.4f  %8.5f   %.4f rad\n',  'W tensión máxima',     x_pol(ip_Wt), y_pol(ip_Wt), L0t_v(ip_Wt), Wt_max)
fprintf('%-28s  %7.4f  %7.4f  %8.5f   %.4f rad\n',  'W compresión máxima',  x_pol(ip_Wc), y_pol(ip_Wc), L0c_v(ip_Wc), Wc_max)
fprintf('%-28s  %7.4f  %7.4f  %8.5f   %.4f rad/m\n','Ef tensión máxima',    x_pol(ip_Eft),y_pol(ip_Eft),L0t_v(ip_Eft),Eft_max)
fprintf('%-28s  %7.4f  %7.4f  %8.5f   %.4f rad/m\n','Ef compresión máxima', x_pol(ip_Efc),y_pol(ip_Efc),L0c_v(ip_Efc),Efc_max)
fprintf('%s\n',repmat('=',1,72))

%% ═══════════════════════════════════════════════════════════
%  GRÁFICAS
%% ═══════════════════════════════════════════════════════════

%% ── Animación ────────────────────────────────────────────────────────
figure
for i = 1:size(P,2)
    clf; hold on

    O2 = [P(1,i);P(2,i)] + Rot(P(3,i))*[-r1/2; 0];
    O4 = [P(1,i);P(2,i)] + Rot(P(3,i))*[ r1/2; 0];
    A  = [P(4,i);P(5,i)] + Rot(P(6,i))*[ r2/2; 0];

    % Centro de barra 3 = punto de conexión con barra 4
    G3 = [P(7,i); P(8,i)];

    % Extremo libre de barra 3
    tip3_i = traj_tip3(:, i);

    % Extremo de barra 4 que conecta con barra 3 (en G3)
    B = G3;   % barra 4 conecta al centro de barra 3

    % Barra 4: de O4 a B (= G3)
    O4_end = [P(10,i);P(11,i)] + Rot(P(12,i))*[ r4/2; 0];   % = B ≈ G3

    Ps2 = traj_s2(:,i);
    Ps4 = traj_s4(:,i);

    %  ── Barras ────────────────────────────────────────────────────────
    line([O2(1) O4(1)], [O2(2) O4(2)], 'Color','k','LineWidth',2)   % barra 1
    line([O2(1)  A(1)], [O2(2)  A(2)], 'Color','r','LineWidth',2)   % barra 2
    line([ A(1) tip3_i(1)], [ A(2) tip3_i(2)], 'Color','g','LineWidth',2)  % barra 3 completa
    line([O4(1)  B(1)], [O4(2)  B(2)], 'Color','b','LineWidth',2)   % barra 4

    % Marca en el punto de conexión B (= centro barra 3)
    plot(B(1), B(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 0.5])

    % Centros de masa
    plot(P(1,i), P(2,i),  'ok','MarkerFaceColor','k','MarkerSize',6)
    plot(P(4,i), P(5,i),  'or','MarkerFaceColor','r','MarkerSize',6)
    plot(P(7,i), P(8,i),  'og','MarkerFaceColor','g','MarkerSize',6)
    plot(P(10,i),P(11,i), 'ob','MarkerFaceColor','b','MarkerSize',6)

    % Trayectoria acumulada del extremo libre de barra 3
    plot(traj_tip3(1,1:i), traj_tip3(2,1:i), 'm-', 'LineWidth', 1.8)
    plot(tip3_i(1), tip3_i(2), 'm^', ...
         'MarkerSize', 10, 'MarkerFaceColor','m', 'MarkerEdgeColor','k')

    % Resorte
    plot(traj_s2(1,1:i), traj_s2(2,1:i), 'r--', 'LineWidth',1.0)
    plot(traj_s4(1,1:i), traj_s4(2,1:i), 'b--', 'LineWidth',1.0)
    plot(Ps2(1),Ps2(2),'s','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r')
    plot(Ps4(1),Ps4(2),'s','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b')
    line([Ps2(1) Ps4(1)],[Ps2(2) Ps4(2)],'Color',[1 0.5 0],'LineWidth',2.5)

    hold off
    axis([-0.15  0.35  -0.15  0.35])   % rango ampliado para ver el extremo
    grid on
    title(sprintf('t = %.2f s  |  k=%d N/m  |  L_0=%.4f m', time(i), k, L0))
    legend('Barra 1','Barra 2','Barra 3','Barra 4','Junta B', ...
           '','','','', ...
           'Tray. ext. B3','Extremo B3', ...
           'Location','northeast','FontSize',7)
    pause(0.05)
end

% %% ── Torque, energía, elongación ──────────────────────────────────────
figure; grid on
plot(P(6,:), LBA(12,:),'r','LineWidth',1.5)
title('Torque vs \phi_2'); xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]')

figure; grid on
plot(time, LBA(12,:),'r','LineWidth',1.5)
title('Torque del motor'); xlabel('Tiempo [s]'); ylabel('Torque [N·m]')

figure; grid on
plot(P(6,:), energia,'b','LineWidth',1.5)
title('Energía potencial del resorte'); xlabel('Tiempo [s]'); ylabel('Energía [J]')

figure; grid on
plot(P(6,:), elongacion,'r','LineWidth',1.5)
title('Elongación del resorte'); xlabel('Tiempo [s]'); ylabel('Elongación [m]')
yline(0,'k--','Longitud natural')

% %% ── Elongación y energía vs L0 ───────────────────────────────────────
figure; hold on; grid on
leg_str=cell(1,N_L0);
for iL=1:N_L0
    plot(phi2_vec,L_vec-L0_sweep(iL),'Color',colores(iL,:),'LineWidth',1.5)
    leg_str{iL}=sprintf('L_0 = %.4f m',L0_sweep(iL));
end
yline(0,'k--','LineWidth',1.2)
xlabel('\phi_2 [rad]'); ylabel('\delta = L - L_0 [m]')
title('Elongación del resorte para cada L_0')
legend(leg_str,'Location','best','FontSize',7); hold off

figure; hold on; grid on
for iL=1:N_L0
    plot(phi2_vec,0.5*k*(L_vec-L0_sweep(iL)).^2,'Color',colores(iL,:),'LineWidth',1.5)
end
xlabel('\phi_2 [rad]'); ylabel('E = \frac{1}{2}k\delta^2 [J]')
title('Energía potencial del resorte para cada L_0')
legend(leg_str,'Location','best','FontSize',7); hold off
% 
% %% ── Rigidez negativa vs L0 (torque) ─────────────────────────────────
fig_tension=figure; hold on; grid on; title('Rigidez negativa TENSIÓN (\delta>0)')
xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]')
fig_compres=figure; hold on; grid on; title('Rigidez negativa COMPRESIÓN (\delta<0)')
xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]')

for iL=1:N_L0
    L0_i=L0_sweep(iL); clr=colores(iL,:);
    TT=TT_store(iL,:);
    dT=diff(TT); Ke=zeros(1,N_pasos-1); Ke(mask_v)=dT(mask_v)./dphi2_v(mask_v);
    dv=L_vec-L0_i;
    mkt=Ke<0 & dv(1:end-1)>0; mkc=Ke<0 & dv(1:end-1)<0;
    figure(fig_tension)
    plot(phi2_vec,TT,'-','Color',[clr 0.3],'LineWidth',0.8)
    for s=1:numel(find(mkt)), ii=find(mkt); ii=ii(s);
        plot(phi2_vec(ii:ii+1),TT(ii:ii+1),'-','Color',clr,'LineWidth',3); end
    figure(fig_compres)
    plot(phi2_vec,TT,'-','Color',[clr 0.3],'LineWidth',0.8)
    for s=1:numel(find(mkc)), ii=find(mkc); ii=ii(s);
        plot(phi2_vec(ii:ii+1),TT(ii:ii+1),'-','Color',clr,'LineWidth',3); end
end
h_t=gobjects(N_L0,1); h_c=gobjects(N_L0,1);
for iL=1:N_L0
    figure(fig_tension); h_t(iL)=plot(nan,nan,'-','Color',colores(iL,:),'LineWidth',2);
    figure(fig_compres); h_c(iL)=plot(nan,nan,'-','Color',colores(iL,:),'LineWidth',2);
end
figure(fig_tension); legend(h_t,leg_str,'Location','best','FontSize',7); hold off
figure(fig_compres); legend(h_c,leg_str,'Location','best','FontSize',7); hold off

figure; hold on; grid on
plot(L0_sweep,W_tension,'r-o','LineWidth',2,'DisplayName','Tensión')
plot(L0_sweep,W_compres,'b-o','LineWidth',2,'DisplayName','Compresión')
plot(L0_sweep,W_vec,'k-o','LineWidth',2,'DisplayName','Total')
xline(L_min,'k--','L_{min}','LabelVerticalAlignment','bottom')
xline(L_max,'k--','L_{max}','LabelVerticalAlignment','bottom')
xlabel('L_0 [m]'); ylabel('Ancho rigidez negativa [\Delta\phi_2 rad]')
title('Ancho rigidez negativa: Tensión vs Compresión vs Total')
legend('Location','best'); hold off
% 
% %% ── Separación horizontal y vertical ─────────────────────────────────% Gráfica: separación horizontal – punto W_comp max
figure_sep_h = figure; hold on; grid on
title(sprintf('Rigidez negativa vs separación horizontal — W_{comp} max, L_0=%.4f m', L0_opt_comp))
xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]')
for is = 1:N_sep
    clr = colores_sep(is,:);
    TT = TT_hor_comp(is,:);
    % Recalcular rigidez efectiva y máscara de K<0 (o usar la ya calculada)
    dT = diff(TT); Ke = zeros(1,N_pasos-1); Ke(mask_v)=dT(mask_v)./dphi2_v(mask_v);
    % Solo regiones de compresión (puedes usar reg_h que ya calculaste)
    % Aquí usamos una máscara simple de K<0 & delta<0 (recalcula delta)
    % (Mejor: almacenaste reg_h, pero si no, calculamos delta)
    % Calculamos elongación para este TT
    % Necesitas recalcular L para este caso... Pero podemos guardar la máscara de región.
    % Para simplificar, podemos usar solo Ke<0 sin filtrar modo.
    nm = Ke<0;
    plot(phi2_vec, TT, '-', 'Color', [clr 0.3], 'LineWidth', 0.8)
    for s = 1:numel(find(nm))
        ii = find(nm); ii=ii(s);
        plot(phi2_vec(ii:ii+1), TT(ii:ii+1), '-', 'Color', clr, 'LineWidth', 3)
    end
end
h=gobjects(N_sep,1);
for is=1:N_sep, h(is)=plot(nan,nan,'-','Color',colores_sep(is,:),'LineWidth',2); end
legend(h, arrayfun(@(d)sprintf('d=%.3f m',d),sep_vec,'UniformOutput',false), 'Location','best','FontSize',7)
hold off

% Repite para vertical, y para el caso de eficiencia...
% 
% %% ── Mapas de calor ───────────────────────────────────────────────────

O20 = [P(1,1);P(2,1)] + Rot(P(3,1))*[-r1/2;0];
O40 = [P(1,1);P(2,1)] + Rot(P(3,1))*[r1/2;0];
A0 = [P(4,1);P(5,1)] + Rot(P(6,1))*[r2/2;0];

G30 = [P(7,1);P(8,1)];
B0 = G30;

O40_end = [P(10,1);P(11,1)] + Rot(P(12,1))*[r4/2/2];



civ_args={'c--','LineWidth',1.5};

figure; imagesc(x_grid,y_grid,mapa_W_tension); set(gca,'YDir','normal')
colorbar; colormap(hot(16384))
title('Ancho rigidez negativa TENSIÓN — L_0 óptimo [\Delta\phi_2 rad]')
xlabel('x [m]'); ylabel('y [m]'); axis([-0.2 0.2 -0.3 0.3]); hold on
plot(CIV24(1,:),CIV24(2,:),civ_args{:})
plot(CIV24(1,1),CIV24(2,1),'c*','MarkerSize',14,'LineWidth',2)
plot(x_grid(ix_t),y_grid(iy_t),'g*','MarkerSize',14,'LineWidth',2)
line([O20(1) O40(1)], [O20(2) O40(2)], 'Color','k','LineWidth',2)   % barra 1
line([O20(1)  A0(1)], [O20(2)  A0(2)], 'Color','r','LineWidth',2)   % barra 2
line([A0(1) tip3_i(1)], [A0(2) tip3_i(2)], 'Color','g','LineWidth',2)  % barra 3 completa
line([O40(1)  B0(1)], [O40(2)  B0(2)], 'Color','b','LineWidth',2)   % barra 4
text(x_grid(ix_t)+.005,y_grid(iy_t)+.01,sprintf('W_{max}=%.4f rad\nL_0^*=%.4f m\n(%.4f,%.4f)',W_t_max,mapa_L0_t(iy_t,ix_t),x_grid(ix_t),y_grid(iy_t)),'Color','g','FontSize',8,'FontWeight','bold')
hold off

figure; imagesc(x_grid,y_grid,mapa_W_compres); set(gca,'YDir','normal')
colorbar; colormap(hot(16384))
title('Ancho rigidez negativa COMPRESIÓN — L_0 óptimo [\Delta\phi_2 rad]')
xlabel('x [m]'); ylabel('y [m]'); axis([-0.6 0.2 -0.3 0.3]); hold on
plot(CIV24(1,:),CIV24(2,:),civ_args{:})
plot(CIV24(1,1),CIV24(2,1),'c*','MarkerSize',14,'LineWidth',2)
plot(x_grid(ix_c),y_grid(iy_c),'g*','MarkerSize',14,'LineWidth',2)
line([O20(1) O40(1)], [O20(2) O40(2)], 'Color','k','LineWidth',2)   % barra 1
line([O20(1)  A0(1)], [O20(2)  A0(2)], 'Color','r','LineWidth',2)   % barra 2
line([A0(1) tip3_i(1)], [A0(2) tip3_i(2)], 'Color','g','LineWidth',2)  % barra 3 completa
line([O40(1)  B0(1)], [O40(2)  B0(2)], 'Color','b','LineWidth',2)   % barra 4
text(x_grid(ix_c)+.005,y_grid(iy_c)+.01,sprintf('W_{max}=%.4f rad\nL_0^*=%.4f m\n(%.4f,%.4f)',W_c_max,mapa_L0_c(iy_c,ix_c),x_grid(ix_c),y_grid(iy_c)),'Color','g','FontSize',8,'FontWeight','bold')
hold off

figure; imagesc(x_grid,y_grid,mapa_efic_t); set(gca,'YDir','normal')
colorbar; colormap(hot(16384))
title('Eficiencia TENSIÓN: W_{tens}/L_{max} [rad/m]')
xlabel('x [m]'); ylabel('y [m]'); axis([-0.2 0.2 -0.3 0.3]); hold on
plot(CIV24(1,:),CIV24(2,:),civ_args{:})
plot(CIV24(1,1),CIV24(2,1),'c*','MarkerSize',14,'LineWidth',2)
plot(x_grid(ix_et),y_grid(iy_et),'g*','MarkerSize',14,'LineWidth',2)
line([O20(1) O40(1)], [O20(2) O40(2)], 'Color','k','LineWidth',2)   % barra 1
line([O20(1)  A0(1)], [O20(2)  A0(2)], 'Color','r','LineWidth',2)   % barra 2
line([A0(1) tip3_i(1)], [A0(2) tip3_i(2)], 'Color','g','LineWidth',2)  % barra 3 completa
line([O40(1)  B0(1)], [O40(2)  B0(2)], 'Color','b','LineWidth',2)   % barra 4
text(x_grid(ix_et)+.005,y_grid(iy_et)+.01,sprintf('Ef_{max}=%.3f rad/m\nL_0^*=%.1f%%\n(%.4f,%.4f)',efic_t_max,mapa_L0pct_t(iy_et,ix_et),x_grid(ix_et),y_grid(iy_et)),'Color','g','FontSize',8,'FontWeight','bold')
hold off

figure; imagesc(x_grid,y_grid,mapa_efic_c); set(gca,'YDir','normal')
colorbar; colormap(hot(16384))
title('Eficiencia COMPRESIÓN: W_{comp}/L_{max} [rad/m]')
xlabel('x [m]'); ylabel('y [m]'); axis([-0.2 0.2 -0.3 0.3]); hold on
plot(CIV24(1,:),CIV24(2,:),civ_args{:})
plot(CIV24(1,1),CIV24(2,1),'c*','MarkerSize',14,'LineWidth',2)
plot(x_grid(ix_ec),y_grid(iy_ec),'g*','MarkerSize',14,'LineWidth',2)
line([O20(1) O40(1)], [O20(2) O40(2)], 'Color','k','LineWidth',2)   % barra 1
line([O20(1)  A0(1)], [O20(2)  A0(2)], 'Color','r','LineWidth',2)   % barra 2
line([A0(1) tip3_i(1)], [A0(2) tip3_i(2)], 'Color','g','LineWidth',2)  % barra 3 completa
line([O40(1)  B0(1)], [O40(2)  B0(2)], 'Color','b','LineWidth',2)   % barra 4
text(x_grid(ix_ec)+.005,y_grid(iy_ec)+.01,sprintf('Ef_{max}=%.3f rad/m\nL_0^*=%.1f%%\n(%.4f,%.4f)',efic_c_max,mapa_L0pct_c(iy_ec,ix_ec),x_grid(ix_ec),y_grid(iy_ec)),'Color','g','FontSize',8,'FontWeight','bold')
hold off
% 
% %% ── Detalle de los 4 casos óptimos ───────────────────────────────────
for ic = 1:4
    W_tot = sum(abs(dphi2_v(region_casos(ic,:))));
    W_t_i = sum(abs(dphi2_v(mask_t_casos(ic,:))));
    W_c_i = sum(abs(dphi2_v(mask_c_casos(ic,:))));
    Lmk   = mapa_Lmax(casos(ic).iy, casos(ic).ix);
    ef_i  = W_tot/max(Lmk,1e-10);
    p_k   = [x_grid(casos(ic).ix); y_grid(casos(ic).iy)];

    figure('Name',sprintf('Caso %d - %s',ic,casos(ic).nombre),'NumberTitle','off')

    subplot(3,1,1); hold on; grid on
    h1=plot(phi2_vec,Torque_base,'b--','LineWidth',1.5);
    h2=plot(phi2_vec,TC(ic,:),'r','LineWidth',2);
    yline(0,'k:','LineWidth',1)
    idx_t2=find(mask_t_casos(ic,:));
    for s=1:numel(idx_t2), ii=idx_t2(s);
        hT=plot(phi2_vec(ii:ii+1),TC(ic,ii:ii+1),'g-','LineWidth',5); end
    idx_c2=find(mask_c_casos(ic,:));
    for s=1:numel(idx_c2), ii=idx_c2(s);
        hC=plot(phi2_vec(ii:ii+1),TC(ic,ii:ii+1),'m-','LineWidth',5); end
    if ~isnan(phi_eq_casos(ic))
        Teq=interp1(phi2_vec,TC(ic,:),phi_eq_casos(ic),'linear','extrap');
        plot(phi_eq_casos(ic),Teq,'k*','MarkerSize',14,'LineWidth',2)
        xline(phi_eq_casos(ic),'k:','LineWidth',1.2)
    end
    title(sprintf('Caso %d: %s\n(%.4f,%.4f) m | L_0=%.4f m\nW=%.4f rad (T:%.4f C:%.4f) Ef=%.3f rad/m', ...
          ic,casos(ic).nombre,p_k(1),p_k(2),casos(ic).L0,W_tot,W_t_i,W_c_i,ef_i))
    lh=[h1,h2]; ls={'Sin resorte','Con resorte'};
    if any(mask_t_casos(ic,:)), lh(end+1)=hT; ls{end+1}='K_{neg} tensión'; end
    if any(mask_c_casos(ic,:)), lh(end+1)=hC; ls{end+1}='K_{neg} compresión'; end
    legend(lh,ls,'Location','best','FontSize',8)
    xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]'); hold off

    subplot(3,1,2); hold on; grid on
    plot(phi2_vec,DK(ic,:),'Color',[0.2 0.6 1],'LineWidth',2)
    for s=1:numel(idx_t2), ii=idx_t2(s);
        plot(phi2_vec(ii:ii+1),DK(ic,ii:ii+1),'g-','LineWidth',5); end
    for s=1:numel(idx_c2), ii=idx_c2(s);
        plot(phi2_vec(ii:ii+1),DK(ic,ii:ii+1),'m-','LineWidth',5); end
    yline(0,'k--','LineWidth',1.2)
    if ~isnan(phi_eq_casos(ic))
        xline(phi_eq_casos(ic),'k:','Equilibrio inestable','LineWidth',1.2,'LabelVerticalAlignment','bottom'); end
    xlabel('\phi_2 [rad]'); ylabel('\delta = L-L_0 [m]')
    title(sprintf('Elongación — L_0=%.4f m',casos(ic).L0)); hold off

    subplot(3,1,3); hold on; grid on
    plot(phi2_vec,EK(ic,:),'Color',[1 0.5 0],'LineWidth',2)
    for s=1:numel(idx_t2), ii=idx_t2(s);
        plot(phi2_vec(ii:ii+1),EK(ic,ii:ii+1),'g-','LineWidth',5); end
    for s=1:numel(idx_c2), ii=idx_c2(s);
        plot(phi2_vec(ii:ii+1),EK(ic,ii:ii+1),'m-','LineWidth',5); end
    if ~isnan(phi_eq_casos(ic))
        Eeq=interp1(phi2_vec,EK(ic,:),phi_eq_casos(ic),'linear','extrap');
        plot(phi_eq_casos(ic),Eeq,'k*','MarkerSize',14,'LineWidth',2)
        xline(phi_eq_casos(ic),'k:','Equilibrio inestable','LineWidth',1.2,'LabelVerticalAlignment','bottom'); end
    xlabel('\phi_2 [rad]'); ylabel('E = \frac{1}{2}k\delta^2 [J]')
    title(sprintf('Energía potencial — L_0=%.4f m',casos(ic).L0)); hold off
    drawnow
end


%% ── GRÁFICAS DE SEPARACIÓN ─────────────────────────────────────────
% Figura 1: Punto de máximo W en compresión
figure('Name','Separación - W_{comp} máximo','NumberTitle','off');
plot(sep_vec*1000, W_hor_comp, 'b-o', 'LineWidth',2, 'DisplayName','Horizontal');
hold on
plot(sep_vec*1000, W_ver_comp, 'r-s', 'LineWidth',2, 'DisplayName','Vertical');
xlabel('Separación d [mm]'); ylabel('Ancho W [rad]');
title(sprintf('W_{comp} max (%.3f, %.3f) m, L_0 = %.4f m', p_opt_comp(1), p_opt_comp(2), L0_opt_comp));
legend('Location','best'); grid on;

% Figura 2: Punto de máxima eficiencia en compresión
figure('Name','Separación - Ef_{comp} máxima','NumberTitle','off');
plot(sep_vec*1000, W_hor_comp_ef, 'b-o', 'LineWidth',2, 'DisplayName','Horizontal');
hold on
plot(sep_vec*1000, W_ver_comp_ef, 'r-s', 'LineWidth',2, 'DisplayName','Vertical');
xlabel('Separación d [mm]'); ylabel('Ancho W [rad]');
title(sprintf('Ef_{comp} max (%.3f, %.3f) m, L_0 = %.4f m', p_opt_comp_ef(1), p_opt_comp_ef(2), L0_opt_comp_ef));
legend('Location','best'); grid on;

%% ── GRÁFICAS TORQUE vs φ₂ PARA LAS SEPARACIONES (4 FIGURAS) ─────────
colores_sep = parula(N_sep);
leg_str = arrayfun(@(d) sprintf('d = %.1f mm', d*1000), sep_vec, 'UniformOutput', false);

% ----- 1. W_comp max – Separación HORIZONTAL -----
figure('Name','Torque vs \phi_2 - W_{comp} max (Horizontal)','NumberTitle','off');
hold on; grid on;
title(sprintf('W_{comp} max (%.3f,%.3f) m, L_0=%.4f m - Horizontal', p_opt_comp(1), p_opt_comp(2), L0_opt_comp));
xlabel('\phi_2 [rad]'); ylabel('Torque motor [N·m]');
plot(phi2_vec, Torque_base, 'k--', 'LineWidth', 1.5);
for is = 1:N_sep
    clr = colores_sep(is,:);
    TT = TT_hor_comp(is,:);
    reg = reg_hor_comp(is,:);
    delta = delta_hor_comp(is,:);
    plot(phi2_vec, TT, '-', 'Color', [clr 0.3], 'LineWidth', 0.8);
    idx_neg = find(reg);
    for s = 1:length(idx_neg)
        ii = idx_neg(s);
        if delta(ii) < 0 || delta(ii+1) < 0
            plot(phi2_vec(ii:ii+1), TT(ii:ii+1), '-', 'Color', clr, 'LineWidth', 3);
        end
    end
end
h_leg = gobjects(N_sep,1);
for is=1:N_sep, h_leg(is)=plot(nan,nan,'-','Color',colores_sep(is,:),'LineWidth',2); end
legend(h_leg, leg_str, 'Location','eastoutside', 'FontSize',7);
hold off

% ----- 2. W_comp max – Separación VERTICAL -----
figure('Name','Torque vs \phi_2 - W_{comp} max (Vertical)','NumberTitle','off');
hold on; grid on;
title(sprintf('W_{comp} max (%.3f,%.3f) m, L_0=%.4f m - Vertical', p_opt_comp(1), p_opt_comp(2), L0_opt_comp));
xlabel('\phi_2 [rad]'); ylabel('Torque motor [N·m]');
plot(phi2_vec, Torque_base, 'k--', 'LineWidth', 1.5);
for is = 1:N_sep
    clr = colores_sep(is,:);
    TT = TT_ver_comp(is,:);
    reg = reg_ver_comp(is,:);
    delta = delta_ver_comp(is,:);
    plot(phi2_vec, TT, '-', 'Color', [clr 0.3], 'LineWidth', 0.8);
    idx_neg = find(reg);
    for s = 1:length(idx_neg)
        ii = idx_neg(s);
        if delta(ii) < 0 || delta(ii+1) < 0
            plot(phi2_vec(ii:ii+1), TT(ii:ii+1), '-', 'Color', clr, 'LineWidth', 3);
        end
    end
end
h_leg = gobjects(N_sep,1);
for is=1:N_sep, h_leg(is)=plot(nan,nan,'-','Color',colores_sep(is,:),'LineWidth',2); end
legend(h_leg, leg_str, 'Location','eastoutside', 'FontSize',7);
hold off

% ----- 3. Ef_comp max (métrica) – Separación HORIZONTAL -----
figure('Name','Torque vs \phi_2 - Ef_{comp} max (Horizontal)','NumberTitle','off');
hold on; grid on;
title(sprintf('Ef_{comp} max (%.3f,%.3f) m, L_0=%.4f m - Horizontal', p_opt_comp_ef(1), p_opt_comp_ef(2), L0_opt_comp_ef));
xlabel('\phi_2 [rad]'); ylabel('Torque motor [N·m]');
plot(phi2_vec, Torque_base, 'k--', 'LineWidth', 1.5);
for is = 1:N_sep
    clr = colores_sep(is,:);
    TT = TT_hor_comp_ef(is,:);
    reg = reg_hor_comp_ef(is,:);
    delta = delta_hor_comp_ef(is,:);
    plot(phi2_vec, TT, '-', 'Color', [clr 0.3], 'LineWidth', 0.8);
    idx_neg = find(reg);
    for s = 1:length(idx_neg)
        ii = idx_neg(s);
        if delta(ii) < 0 || delta(ii+1) < 0
            plot(phi2_vec(ii:ii+1), TT(ii:ii+1), '-', 'Color', clr, 'LineWidth', 3);
        end
    end
end
h_leg = gobjects(N_sep,1);
for is=1:N_sep, h_leg(is)=plot(nan,nan,'-','Color',colores_sep(is,:),'LineWidth',2); end
legend(h_leg, leg_str, 'Location','eastoutside', 'FontSize',7);
hold off

% ----- 4. Ef_comp max (métrica) – Separación VERTICAL -----
figure('Name','Torque vs \phi_2 - Ef_{comp} max (Vertical)','NumberTitle','off');
hold on; grid on;
title(sprintf('Ef_{comp} max (%.3f,%.3f) m, L_0=%.4f m - Vertical', p_opt_comp_ef(1), p_opt_comp_ef(2), L0_opt_comp_ef));
xlabel('\phi_2 [rad]'); ylabel('Torque motor [N·m]');
plot(phi2_vec, Torque_base, 'k--', 'LineWidth', 1.5);
for is = 1:N_sep
    clr = colores_sep(is,:);
    TT = TT_ver_comp_ef(is,:);
    reg = reg_ver_comp_ef(is,:);
    delta = delta_ver_comp_ef(is,:);
    plot(phi2_vec, TT, '-', 'Color', [clr 0.3], 'LineWidth', 0.8);
    idx_neg = find(reg);
    for s = 1:length(idx_neg)
        ii = idx_neg(s);
        if delta(ii) < 0 || delta(ii+1) < 0
            plot(phi2_vec(ii:ii+1), TT(ii:ii+1), '-', 'Color', clr, 'LineWidth', 3);
        end
    end
end
h_leg = gobjects(N_sep,1);
for is=1:N_sep, h_leg(is)=plot(nan,nan,'-','Color',colores_sep(is,:),'LineWidth',2); end
legend(h_leg, leg_str, 'Location','eastoutside', 'FontSize',7);
hold off


% 
% %% ── Poloida coloreada por W_best ─────────────────────────────────────
figure('Name','Poloida — W mejor (sin métrica)','NumberTitle','off','Position',[100 100 900 750])
subplot(2,1,1); hold on; grid on; axis equal
scatter(x_pol,y_pol,40,Wt_v,'filled'); colormap(hot(16384))
cb=colorbar; cb.Label.String='W_{tensión} [rad]'; clim([0 max(Wt_v)+1e-10])
plot(x_pol(ip_Wt),y_pol(ip_Wt),'c*','MarkerSize',16,'LineWidth',2)
text(x_pol(ip_Wt)+.004,y_pol(ip_Wt)+.005,sprintf('W_{max}=%.4f rad\nL_0^*=%.4f m\n(%.4f,%.4f)',Wt_max,L0t_v(ip_Wt),x_pol(ip_Wt),y_pol(ip_Wt)),'Color','c','FontSize',8,'FontWeight','bold')
xlabel('x [m]'); ylabel('y [m]'); title('Poloida CIV 2-4 — Mejor W_{tensión} [\Delta\phi_2 rad]'); hold off

subplot(2,1,2); hold on; grid on; axis equal
scatter(x_pol,y_pol,40,Wc_v,'filled'); colormap(hot(16384))
cb=colorbar; cb.Label.String='W_{compresión} [rad]'; clim([0 max(Wc_v)+1e-10])
plot(x_pol(ip_Wc),y_pol(ip_Wc),'c*','MarkerSize',16,'LineWidth',2)
text(x_pol(ip_Wc)+.004,y_pol(ip_Wc)+.005,sprintf('W_{max}=%.4f rad\nL_0^*=%.4f m\n(%.4f,%.4f)',Wc_max,L0c_v(ip_Wc),x_pol(ip_Wc),y_pol(ip_Wc)),'Color','c','FontSize',8,'FontWeight','bold')
xlabel('x [m]'); ylabel('y [m]'); title('Poloida CIV 2-4 — Mejor W_{compresión} [\Delta\phi_2 rad]'); hold off

%% ── Poloida coloreada por Ef_best ────────────────────────────────────
figure('Name','Poloida — Eficiencia (con métrica)','NumberTitle','off','Position',[150 150 900 750])
subplot(2,1,1); hold on; grid on; axis equal
scatter(x_pol,y_pol,40,Eft_v,'filled'); colormap(hot(16384))
cb=colorbar; cb.Label.String='Ef_{tensión} [rad/m]'; clim([0 max(Eft_v)+1e-10])
plot(x_pol(ip_Eft),y_pol(ip_Eft),'c*','MarkerSize',16,'LineWidth',2)
text(x_pol(ip_Eft)+.004,y_pol(ip_Eft)+.005,sprintf('Ef_{max}=%.3f rad/m\nL_0^*=%.4f m\n(%.4f,%.4f)',Eft_max,L0t_v(ip_Eft),x_pol(ip_Eft),y_pol(ip_Eft)),'Color','c','FontSize',8,'FontWeight','bold')
xlabel('x [m]'); ylabel('y [m]'); title('Poloida CIV 2-4 — Mejor Ef_{tensión} = W/L_{max} [rad/m]'); hold off

subplot(2,1,2); hold on; grid on; axis equal
scatter(x_pol,y_pol,40,Efc_v,'filled'); colormap(hot(16384))
cb=colorbar; cb.Label.String='Ef_{compresión} [rad/m]'; clim([0 max(Efc_v)+1e-10])
plot(x_pol(ip_Efc),y_pol(ip_Efc),'c*','MarkerSize',16,'LineWidth',2)
text(x_pol(ip_Efc)+.004,y_pol(ip_Efc)+.005,sprintf('Ef_{max}=%.3f rad/m\nL_0^*=%.4f m\n(%.4f,%.4f)',Efc_max,L0c_v(ip_Efc),x_pol(ip_Efc),y_pol(ip_Efc)),'Color','c','FontSize',8,'FontWeight','bold')
xlabel('x [m]'); ylabel('y [m]'); title('Poloida CIV 2-4 — Mejor Ef_{compresión} = W/L_{max} [rad/m]'); hold off

%% ═══════════════════════════════════════════════════════════
%  FUNCIONES LOCALES
%% ═══════════════════════════════════════════════════════════
function region = calcular_region_eq(Torque, K_eff, neg_mask, L_vec, L0, k, N_pasos, modo)
%CALCULAR_REGION_EQ  Región contigua de K_eff<0 alrededor del equilibrio
%  inestable del modo indicado ('tension' o 'compresion').
    region  = false(1, N_pasos-1);
    sign_ch = find(diff(sign(Torque)));
    if isempty(sign_ch), return; end
 
    E_v  = 0.5*k*(L_vec - L0).^2;
    dv   = L_vec - L0;
    n_sc = numel(sign_ch);
    E_cross    = zeros(1, n_sc);
    Kneg_cross = false(1, n_sc);
 
    % PASO 1 — calcular E_cross y Kneg_cross
    for sc = 1:n_sc
        idx_sc = sign_ch(sc);
        E_cross(sc) = (E_v(idx_sc) + E_v(min(idx_sc+1, N_pasos))) / 2;
        if idx_sc <= N_pasos-1
            Kneg_cross(sc) = neg_mask(idx_sc);
        end
    end
 
    % PASO 2 — copiar a E_cand
    E_cand = E_cross;
 
       % PASO 4 primero — eliminar cruces del modo contrario
    for sc = 1:n_sc
        d_cross = (dv(sign_ch(sc)) + dv(min(sign_ch(sc)+1, N_pasos))) / 2;
        if strcmp(modo,'tension')    && d_cross <  0, E_cand(sc) = -Inf; end
        if strcmp(modo,'compresion') && d_cross >= 0, E_cand(sc) = -Inf; end
    end
    
    if all(E_cand == -Inf), return; end   % no hay cruces de este modo
    
    % PASO 3 después — entre los válidos del modo, preferir K_eff<0
    mode_y_Kneg = (E_cand > -Inf) & Kneg_cross;
    if any(mode_y_Kneg)
        E_cand(~Kneg_cross) = -Inf;   % solo eliminar si hay K_eff<0 en este modo
    end
 
    % Si todos quedaron en -Inf no hay equilibrio de ese modo
    if all(E_cand == -Inf), return; end
 
    % PASO 5 — elegir el cruce con mayor energía
    [~, bc] = max(E_cand);
    anchor  = min(max(sign_ch(bc), 1), N_pasos-1);
 
    if ~neg_mask(anchor)
        neg_idx = find(neg_mask);
        if isempty(neg_idx), return; end
        [~, ci] = min(abs(neg_idx - anchor));
        anchor  = neg_idx(ci);
    end
 
    % PASO 6 — expandir la región contigua de K_eff < 0
    sr = anchor; er = anchor;
    while sr > 1         && neg_mask(sr-1), sr = sr-1; end
    while er < N_pasos-1 && neg_mask(er+1), er = er+1; end
    region(sr:er) = true;
end

function [Fq, J] = ecuaciones(q, t, r1, r2, r3, r4, phi20, omega20, alpha20)
    Fq = [q(1)-cos(q(3))*r1/2 - q(4)+cos(q(6))*r2/2;   % nodo O2: barra1-barra2
          q(2)-sin(q(3))*r1/2 - q(5)+sin(q(6))*r2/2;
          q(4)+cos(q(6))*r2/2 - q(7)+cos(q(9))*r3/2;    % nodo A: barra2-barra3
          q(5)+sin(q(6))*r2/2 - q(8)+sin(q(9))*r3/2;
          q(7) - q(10) - cos(q(12))*r4/2;                % nodo B: barra3-barra4
          q(8) - q(11) - sin(q(12))*r4/2;                %   barra4 en CENTRO de barra3
          q(1)+cos(q(3))*r1/2 - q(10)+cos(q(12))*r4/2;  % nodo O4: barra1-barra4
          q(2)+sin(q(3))*r1/2 - q(11)+sin(q(12))*r4/2;
          q(1) - r1/2;   % barra 1 fija
          q(2);
          q(3);
          q(6) - 0.5*alpha20*t^2 - omega20*t - phi20];   % ley de movimiento
 
    J = [1  0   sin(q(3))*r1/2   -1  0  -sin(q(6))*r2/2   0  0   0                0  0   0              ;
         0  1  -cos(q(3))*r1/2    0 -1   cos(q(6))*r2/2   0  0   0                0  0   0              ;
         0  0   0                 1  0  -sin(q(6))*r2/2   -1  0  -sin(q(9))*r3/2   0  0   0              ;
         0  0   0                 0  1   cos(q(6))*r2/2    0 -1   cos(q(9))*r3/2   0  0   0              ;
         0  0   0                 0  0   0                  1  0   0               -1  0   sin(q(12))*r4/2;  % ← cambiado
         0  0   0                 0  0   0                  0  1   0                0 -1  -cos(q(12))*r4/2;  % ← cambiado
         1  0  -sin(q(3))*r1/2    0  0   0                  0  0   0               -1  0  -sin(q(12))*r4/2;
         0  1   cos(q(3))*r1/2    0  0   0                  0  0   0                0 -1   cos(q(12))*r4/2;
         1  0   0  0  0  0  0  0  0  0  0  0;
         0  1   0  0  0  0  0  0  0  0  0  0;
         0  0   1  0  0  0  0  0  0  0  0  0;
         0  0   0  0  0  1  0  0  0  0  0  0];
end



function Jp = jacobian_punto(q, v, r1, r2, r3, r4)
Jp = [0 0  v(3)*cos(q(3))*r1/2   0 0 -v(6)*cos(q(6))*r2/2   0 0  0                    0 0  0                  ;
    0 0  v(3)*sin(q(3))*r1/2   0 0 -v(6)*sin(q(6))*r2/2   0 0  0                    0 0  0                  ;
    0 0  0                      0 0 -v(6)*cos(q(6))*r2/2   0 0 -v(9)*cos(q(9))*r3/2  0 0  0                  ;
    0 0  0                      0 0 -v(6)*sin(q(6))*r2/2   0 0 -v(9)*sin(q(9))*r3/2  0 0  0                  ;
    0 0  0   0 0 0   0 0 -v(9)*cos(q(9))*r3/2  0 0  v(12)*cos(q(12))*r4/2;
    0 0  0   0 0 0   0 0 -v(9)*sin(q(9))*r3/2  0 0  v(12)*sin(q(12))*r4/2;
    0 0 -v(3)*cos(q(3))*r1/2   0 0  0   0 0  0   0 0 -v(12)*cos(q(12))*r4/2;
    0 0 -v(3)*sin(q(3))*r1/2   0 0  0   0 0  0   0 0 -v(12)*sin(q(12))*r4/2;
    zeros(4,12)];
end


%% ═══════════════════════════════════════════════════════════
%  FUNCIÓN AUXILIAR: Ancho de rigidez negativa para un par
%  de anclajes inicial y una L0 dada.
%% ═══════════════════════════════════════════════════════════
function W = calcular_W_para_punto(p2_inic, p4_inic, L0, modo)
% Calcula el ancho total de rigidez negativa (en rad de phi2)
% para un resorte cuyos anclajes en t=0 son p2_inic (en barra 2) y
% p4_inic (en barra 4), con longitud natural L0, en el modo dado.
% Usa las variables globales del script principal (P, Torque_base, etc.)
% ── Obtener variables globales necesarias ─────────────────
%global G2_foto G4_foto q_foto P Torque_base J_hist time N_pasos ...
 %   dphi2_v mask_v phi2_vec

% Coordenadas locales constantes
r2_loc = Rot(q_foto(6))'  * (p2_inic - G2_foto);
r4_loc = Rot(q_foto(12))' * (p4_inic - G4_foto);

% Contribución del resorte a lo largo del tiempo
dt = zeros(1, N_pasos);
for i = 1:N_pasos
    ps2 = [P(4,i); P(5,i)] + Rot(P(6,i))  * r2_loc;
    ps4 = [P(10,i); P(11,i)] + Rot(P(12,i)) * r4_loc;
    L = norm(ps4 - ps2);
    if L > 1e-10
        u = (ps4 - ps2) / L;
    else
        u = [0; 0];
    end
    F = 2 * (L - L0) * u;   % k = 2 N/m (¡ajusta si cambiaste k!)
    dR2 = [-sin(P(6,i))  -cos(P(6,i));  cos(P(6,i))  -sin(P(6,i))];
    dR4 = [-sin(P(12,i)) -cos(P(12,i)); cos(P(12,i)) -sin(P(12,i))];
    Q = zeros(12,1);
    Q(4)  =  F(1);
    Q(5)  =  F(2);
    Q(6)  = (dR2 * r2_loc)' * F;
    Q(10) = -F(1);
    Q(11) = -F(2);
    Q(12) = -(dR4 * r4_loc)' * F;
    lambda_i = J_hist(:,:,i)' \ Q;
    dt(i) = -lambda_i(12);   % contribución al par motor
end

T_total = Torque_base + dt;

% Rigidez efectiva
dT = diff(T_total);
K_eff = zeros(1, N_pasos-1);
K_eff(mask_v) = dT(mask_v) ./ dphi2_v(mask_v);
neg_mask = K_eff < 0;

% Elongación
L_vec = zeros(1, N_pasos);
for i = 1:N_pasos
    ps2 = [P(4,i); P(5,i)] + Rot(P(6,i))  * r2_loc;
    ps4 = [P(10,i); P(11,i)] + Rot(P(12,i)) * r4_loc;
    L_vec(i) = norm(ps4 - ps2);
end
delta_vec = L_vec - L0;

% Seleccionar región según el modo
region = calcular_region_eq(T_total, K_eff, neg_mask, L_vec, L0, 2, N_pasos, modo);

% Sumar ancho angular
W = sum(abs(dphi2_v(region & (delta_vec(1:end-1) >= 0))));   % tensión
if strcmp(modo, 'compresion')
    W = sum(abs(dphi2_v(region & (delta_vec(1:end-1) < 0))));
end
end

function m = Rot(theta)
    m = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end














