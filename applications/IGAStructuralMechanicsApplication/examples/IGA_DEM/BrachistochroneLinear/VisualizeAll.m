faces = importdata('faces.txt')
[s,m] = size(faces)
figure;
hold on;
%  embedded_faces = importdata('facesembedded.txt')
%  [s,m] = size(embedded_faces)
% %for n = 1:s
plot3(faces(:,1), faces(:,2), faces(:,3), '.')
% plot3(embedded_faces(:,1), embedded_faces(:,2), embedded_faces(:,3), '.', 'color', 'red')
%end

edges = importdata('edges.txt')
[s,m] = size(edges)
figure;
%for n = 1:s
    plot3(edges(:,1), edges(:,2), edges(:,3), '.')
%end



cedges = importdata('c_edges.txt');
[s,m] = size(cedges);
figure;
hold on;
for n = 1:s
    X = [cedges(n,1) cedges(n,4)];
    Y = [cedges(n,2) cedges(n,5)];
    Z = [cedges(n,3) cedges(n,6)];
    plot3(X, Y, Z, 'black')
end

figure;
hold on;
plot3(faces(:,1), faces(:,2), faces(:,3), '.')
plot3(edges(:,1), edges(:,2), edges(:,3),'.','color','red')
plot3(embedded_faces(:,1), embedded_faces(:,2), embedded_faces(:,3), '.', 'color', 'black')
plot3(cedges(:,1), cedges(:,2), cedges(:,3),'.', 'color','green')
plot3(cedges(:,4), cedges(:,5), cedges(:,6),'.', 'color','black')
for n = 1:s
    X = [cedges(n,1) cedges(n,4)];
    Y = [cedges(n,2) cedges(n,5)];
    Z = [cedges(n,3) cedges(n,6)];
    plot3(X, Y, Z, 'black')
end