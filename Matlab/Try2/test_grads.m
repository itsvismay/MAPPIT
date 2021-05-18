q_i = [0,
        0,
        0,
    1.125,
    1.125,
        1,
     2.25,
     2.25,
        2,
    3.375,
    3.375,
        3,
     5.25,
     3.75,
        4,
      7.5,
     3.75,
        5,
     9.75,
     3.75,
        6,
   11.625,
    3.375,
        7,
    12.75,
     2.25,
        8,
   13.875,
    1.125,
        9,
       15,
        0,
       10]; %3*nodes
   
dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
pv = 1;
K = 1;
e = K*0.5*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv) ).^2;

%         e = e + K*0.5*(q_i(end).^2 - (sum((dx(:,1).^2  + dx(:,2).^2))/(pv*pv)) ).^2;

gradE = ones(size(q_i));
        gradE = gradE*K*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv));
        dedx2 = -(dx(:,1).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;
        dedx1 = (dx(:,1).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;
        dedy2 = -(dx(:,2).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;
        dedy1 = (dx(:,2).*((dx(:,1).^2  + dx(:,2).^2).^(-1/2)))/pv;

        dEdq_left = zeros(numel(q_i)/3, 3);
        dEdq_right = zeros(numel(q_i)/3, 3);
        dEdq_left(1:end-1,1) = dedx1;
        dEdq_left(1:end-1,2) = dedy1;
        dEdq_right(2:end, 1) = dedx2;
        dEdq_right(2:end, 2) = dedy2;

        dEdq = dEdq_left + dEdq_right;
        flatdEdq = reshape(dEdq', size(dEdq,1)*size(dEdq,2), 1);
        gradE(1:end-1) = gradE(1:end-1).*flatdEdq(1:end-1);

dEdq = zeros(size(q_i));
dEdq1 = dEdq;


for i=1:size(q_i,1)
    q_i(i) =q_i(i)+ 1e-4;
    dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
    el = K*0.5*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv) ).^2;
    q_i(i) =q_i(i)- 1e-4;
    q_i(i) =q_i(i)- 1e-4;
    dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
    er = K*0.5*(q_i(end) - (sum(sqrt(dx(:,1).^2  + dx(:,2).^2))/pv) ).^2;
    q_i(i) =q_i(i)+ 1e-4;
    
    dEdq(i) = (el - er)/(2*1e-4);
    %dEdq1(i) = (el - e)/(1e-4);
end
dEdq - gradE

