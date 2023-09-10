function M = estimateCameraProjectionMatrix(impoints, objpoints)
    [point,w] = size(objpoints);
    p = [];
    for i = 1:point
        Xo = objpoints(i,1);
        Yo = objpoints(i,2);
        Zo = objpoints(i,3);

        x = impoints(i,1);
        y = impoints(i,2);

        eq = [-Xo,-Yo,-Zo,-1,0,0,0,0,x*Xo,x*Yo,x*Zo,x;
              0,0,0,0,-Xo,-Yo,-Zo,-1,y*Xo,y*Yo,y*Zo,y];

        p = [p;eq];
    end

    [U,S,V] = svd(p,'econ');
    q = V(:,end);
    M = reshape(q,4,3)';
end