clc; clear all; close all
%% Hooken_torsion_spring_CIV




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

k  = 20;       % N/m
L0 = 0.1480;  % m




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




ktO2 = 0.1;   % Nm/rad
KtA  = 0.05;
ktB  = 0.1;
kt4  = 0.05;


theta_eq_O2 = q_foto(6);
theta_eq_A  = q_foto(9) - q_foto(6);
theta_eq_B  = q_foto(9)  - q_foto(12) ;
theta_eq_O4 = q_foto(12);

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


% Prealocación para resortes de torsión
deform_O2  = zeros(1, ceil(36/0.1)+1);   % deformación angular barra 2 - tierra
deform_A   = zeros(1, ceil(36/0.1)+1);   % deformación en A (θ3-θ2)
deform_B   = zeros(1, ceil(36/0.1)+1);   % deformación en B (θ4-θ3)
deform_O4  = zeros(1, ceil(36/0.1)+1);   % deformación barra 4 - tierra

ener_O2 = zeros(1, ceil(36/0.1)+1);
ener_A  = zeros(1, ceil(36/0.1)+1);
ener_B  = zeros(1, ceil(36/0.1)+1);
ener_O4 = zeros(1, ceil(36/0.1)+1);
ener_tors_total = zeros(1, ceil(36/0.1)+1);

Torque_torsion = zeros(1, ceil(36/0.1)+1);
Torque_lineal  = zeros(1, ceil(36/0.1)+1);
 
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

    
    % --- Resortes de torsión (CORREGIDO) ---
    Q_tors = zeros(12,1);

    % 1. Resorte en O2 (barra 2 – tierra)
    tau_O2 = ktO2 * (theta_eq_O2 - q(6));    % restaurador
    Q_tors(6) = Q_tors(6) + tau_O2;

    % 2. Resorte en A (barra 3 – barra 2)
    delta_theta_A = (q(9) - q(6)) - theta_eq_A;
    tau_A = KtA * delta_theta_A;             % par que actúa sobre barra 2 para reducir la diferencia
    Q_tors(6) = Q_tors(6) + tau_A;           % contribución a θ2: +∂U/∂θ2 = +KtA*(…)
    Q_tors(9) = Q_tors(9) - tau_A;           % contribución a θ3: -∂U/∂θ3 = -KtA*(…)

    % 3. Resorte en B (barra 4 – barra 3)
    delta_theta_B = (q(9) - q(12)) - theta_eq_B;
    tau_B = ktB * delta_theta_B;
    Q_tors(9) = Q_tors(9) + tau_B;           % +∂U/∂θ3
    Q_tors(12) = Q_tors(12) - tau_B;         % -∂U/∂θ4

    % 4. Resorte en O4 (barra 4 – tierra)
    tau_O4 = kt4 * (theta_eq_O4 - q(12));
    Q_tors(12) = Q_tors(12) + tau_O4;


 
    % Lambda total (lineal + torsión)
    lambda           = inv(J') * (Mmat*a - Fext - Q_spring - Q_tors);
    LBA(:,contador)  = lambda;

    % Lambda base (sin ningún resorte)
    lambda_base           = inv(J') * (Mmat*a - Fext);
    Torque_base(contador) = lambda_base(12);

    % Lambda solo con torsión; rigidez positiva
    lambda_torsion        = inv(J') * (Mmat*a - Fext - Q_tors);
    Torque_torsion(contador) = lambda_torsion(12);


    %lambda solo con resorte lineal; rigidez negativa

    lambda_lineal = inv(J') * (Mmat*a - Fext - Q_spring);
    Torque_lineal(contador) = lambda_lineal(12);
    lambda_base           = inv(J') * (Mmat*a - Fext);
    Torque_base(contador) = lambda_base(12);
 
    elongacion(contador) = delta;
    energia(contador)    = 0.5 * k * delta^2;
    longitud(contador)   = L_actual;
    traj_s2(:,contador)  = Ps2;
    traj_s4(:,contador)  = Ps4;

    % Deformaciones y energías de los resortes de torsión
    deform_O2(contador) = q(6) - theta_eq_O2;
    ener_O2(contador)   = 0.5 * ktO2 * deform_O2(contador)^2;

    deform_A(contador) = (q(9) - q(6)) - theta_eq_A;
    ener_A(contador)   = 0.5 * KtA * deform_A(contador)^2;

    deform_B(contador) = (q(9) - q(12)) - theta_eq_B;
    ener_B(contador)   = 0.5 * ktB * deform_B(contador)^2;

    deform_O4(contador) = q(12) - theta_eq_O4;
    ener_O4(contador)   = 0.5 * kt4 * deform_O4(contador)^2;

    ener_tors_total(contador) = ener_O2(contador) + ener_A(contador) + ener_B(contador) + ener_O4(contador);

 
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


save("Hookens_torsion_spring_CIV_25.mat")
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

% Graficas de torque total, torque con linear spring y torque con torsion
% spring

figure; 
plot(P(6,:), LBA(12,:),'r','LineWidth',1.5,'DisplayName','Torque total')
grid on
hold on
plot(P(6,:),Torque_torsion,'g','LineWidth',1.5,'DisplayName','Torque torsion spring')
plot(P(6,:),Torque_lineal,'b','LineWidth',1.5,'DisplayName','Torque linear spring')
legend('Torque total','Torque torsion spring','Torque linear springs')
title('Torque vs \phi_2'); xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]')
hold off

figure;
subplot(3,1,1)
plot(P(6,:), LBA(12,:), 'r'); grid on; title('Torque total')
subplot(3,1,2)
plot(P(6,:), Torque_lineal, 'b'); grid on; title('Torque con resorte lineal')
subplot(3,1,3)
plot(P(6,:), Torque_torsion, 'g'); grid on; title('Torque con resortes de torsión')
xlabel('\phi_2 [rad]')


figure; hold on; grid on
h1 = plot(P(6,:), LBA(12,:), 'r-', 'LineWidth', 2);
h2 = plot(P(6,:), Torque_lineal, 'b--', 'LineWidth', 1.5);
h3 = plot(P(6,:), Torque_torsion, 'g:', 'LineWidth', 2);
legend([h1 h2 h3], {'Total','Solo lineal','Solo torsión'})
xlabel('\phi_2 [rad]'); ylabel('Torque [N·m]')
title('Comparación de pares motores')



% Deformaciones angulares y elongacion

figure; 
plot(P(6,:), elongacion,'r','LineWidth',1.5)
grid on
title('Elongación del resorte'); xlabel('\phi_2[rad]'); ylabel('Elongación [m]')
yline(0,'k--','Longitud natural')

figure;
subplot(2,2,1); plot(P(6,:), deform_O2); grid on; title('Resorte O2'); ylabel('\Delta\theta [rad]');
subplot(2,2,2); plot(P(6,:), deform_A);  grid on; title('Resorte A'); ylabel('\Delta\theta [rad]');
subplot(2,2,3); plot(P(6,:), deform_B);  grid on; title('Resorte B'); ylabel('\Delta\theta [rad]');
subplot(2,2,4); plot(P(6,:), deform_O4); grid on; title('Resorte O4'); ylabel('\Delta\theta [rad]');
sgtitle('Deformación angular – Resortes de torsión');

% Energías
figure; 
plot(P(6,:), energia,'b','LineWidth',1.5)
grid on
title('Energía potencial del resorte'); xlabel('\phi_2[rad]'); ylabel('Energía [J]')

figure;
subplot(2,2,1); plot(P(6,:), ener_O2); grid on; title('Resorte O2'); ylabel('Energía [J]');
subplot(2,2,2); plot(P(6,:), ener_A);  grid on; title('Resorte A'); ylabel('Energía [J]');
subplot(2,2,3); plot(P(6,:), ener_B);  grid on; title('Resorte B'); ylabel('Energía [J]');
subplot(2,2,4); plot(P(6,:), ener_O4); grid on; title('Resorte O4'); ylabel('Energía [J]');
sgtitle('Energía potencial – Resortes de torsión');

%% ── Energía total del sistema ────────────────────────────────────────
figure;
plot(P(6,:), energia, 'b', 'LineWidth', 1.5, 'DisplayName', 'Resorte lineal');
hold on; grid on;
plot(P(6,:), ener_tors_total, 'r', 'LineWidth', 1.5, 'DisplayName', 'Resortes de torsión (suma)');
plot(P(6,:), energia + ener_tors_total, 'k', 'LineWidth', 2, 'DisplayName', 'Energía total');
xlabel('\phi_2 [rad]'); ylabel('Energía [J]');
title('Energía potencial del sistema');
legend('Location', 'best');
hold off;





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
