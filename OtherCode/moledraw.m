function moledraw(xyz_tb,adj_list)
nhc=size(keys(adj_list),1);
for j=1:nhc
    nb_node=adj_list({num2str(j-1)});
    nb_node=sort(nb_node{:});
    nsize=size(nb_node,2);
    for i=1:nsize
        line_parent_child(xyz_tb,j,nb_node(1,i)+1);
        hold on;
    end
end
for i=1:size(xyz_tb,1)
    text(xyz_tb(i,1),xyz_tb(i,2),xyz_tb(i,3),num2str(i));
end
axis equal
grid on
box off
xlabel("X")
ylabel("Y")
zlabel("Z")

end