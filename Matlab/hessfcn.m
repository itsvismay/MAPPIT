function [H] = hessfcn(q_i,lambda)
    global num_agents scene;
    Q = reshape(q_i, numel(q_i)/num_agents, num_agents);    
    [e, H_ke] = kinetic_energy_with_hessian(Q, scene, 1, 1);
    H = zeros(size(q_i,1),size(q_i,1));
    for ii=1:num_agents
        H(1+(ii-1)*size(Q,1):ii*size(Q,1),1+(ii-1)*size(Q,1):ii*size(Q,1)) = H_ke(:,:,ii);
        %H_full(end-ii+1,end-ii+1) = 1; %assign the last block to identity
    end
    
end

function [e, HT] = kinetic_energy_with_hessian(Q, scene, K, constraint_set)
    if sum(K)==0
        e=0;
        HT = zeros(size(Q,1), size(Q,1), size(Q,2));
        return;
    end
    HT = zeros(size(Q,1), size(Q,1), size(Q,2));
    e=0;
    for i=1:numel(scene.agents)

        q_i = Q(:, i); %3*nodes
        m = scene.agents(i).mass;
        dx = reshape(q_i(4:end) - q_i(1:end -3), 3, numel(q_i)/3-1)';
        
        d2Edq2 = zeros(numel(q_i),numel(q_i));
         
        if constraint_set == 1
            %spring energy
            e = e + K*sum(0.5*m*sum(dx(:, 1:3).*dx(:,1:3),2)); %kinetic energy
            
            for kk=1:numel(q_i)/3-1
                d2Edq2((kk-1)*3+1,(kk-1)*3+1) = d2Edq2((kk-1)*3+1,(kk-1)*3+1) + m;
                d2Edq2((kk-1)*3+2,(kk-1)*3+2) = d2Edq2((kk-1)*3+2,(kk-1)*3+2) + m;
                d2Edq2((kk-1)*3+3,(kk-1)*3+3) = d2Edq2((kk-1)*3+3,(kk-1)*3+3) + m;
                
                d2Edq2((kk-1)*3+1,(kk-1)*3+4) = d2Edq2((kk-1)*3+1,(kk-1)*3+4) - m;
                d2Edq2((kk-1)*3+2,(kk-1)*3+5) = d2Edq2((kk-1)*3+2,(kk-1)*3+5) - m;
                d2Edq2((kk-1)*3+3,(kk-1)*3+6) = d2Edq2((kk-1)*3+3,(kk-1)*3+6) - m;
                
                d2Edq2((kk-1)*3+4,(kk-1)*3+1) = d2Edq2((kk-1)*3+4,(kk-1)*3+1) - m;
                d2Edq2((kk-1)*3+5,(kk-1)*3+2) = d2Edq2((kk-1)*3+5,(kk-1)*3+2) - m;
                d2Edq2((kk-1)*3+6,(kk-1)*3+3) = d2Edq2((kk-1)*3+6,(kk-1)*3+3) - m;
                
                d2Edq2((kk-1)*3+4,(kk-1)*3+4) = d2Edq2((kk-1)*3+4,(kk-1)*3+4) + m;
                d2Edq2((kk-1)*3+5,(kk-1)*3+5) = d2Edq2((kk-1)*3+5,(kk-1)*3+5) + m;
                d2Edq2((kk-1)*3+6,(kk-1)*3+6) = d2Edq2((kk-1)*3+6,(kk-1)*3+6) + m;
            end
        end
               
        HT(:,:,i) = d2Edq2;
            
    end
end
