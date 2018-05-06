function[normalized_normal_vectors]=STL_Export(nodes_file, triangles_file, STL_file_name,solid)
%This is a tool to export 3D graphics from a Tri_Surface file to an ASCII STL file.

%Before the export, the Tri_Surface file can be used in Matlab for
%visualization, some calculations or modifications. Thereafter, the data
%for further processing often re-used in the STL format.

%INPUT:

% nodes_file: Is the matrix of 3D - points
% triangles_file: Is the matrix of the triangles
% STL_file_name: Is the name for the generating STL file
% solid: Is the name of the file inside the STL file

%OUTPUT:

% nomalized_normal_vectors: Is the matrix of the normal vectors of each
%triangle.

% And of course: The STL file. 

%Starting the function

%calculate the size of the surface data

tic
node_num=length(nodes_file);
triangle_num=length(triangles_file);

%print the size of the surface data

fprintf (1,'\n');
fprintf (1,'  TRI_SURFACE data:\n');
fprintf (1,'  Number of nodes     = %d\n',node_num);
fprintf (1,'  Number of triangles = %d\n',triangle_num);

%compute the normal vectors without a loop

points_triangles=[nodes_file(triangles_file,1),nodes_file(triangles_file,2),nodes_file(triangles_file,3)];
points_one=points_triangles(1:length(points_triangles)/3,:);
points_two=points_triangles(length(points_triangles)/3+1:length(points_triangles)/3*2,:);
points_three=points_triangles(length(points_triangles)/3*2+1:length(points_triangles),:);
vectors_one=points_two-points_one;
vectors_two=points_three-points_one;
normal_vectors=cross(vectors_one,vectors_two);
norms=repmat(sqrt(sum(normal_vectors.^2,2)),[1,3]);
normalized_normal_vectors=normal_vectors./norms;

%create the output matrix

output=zeros(length(points_one)*4,3);
for i=1:length(points_one)
    output(i*4-3,:)=normalized_normal_vectors(i,:);
    output(i*4-2,:)=points_one(i,:);
    output(i*4-1,:)=points_two(i,:);
    output(i*4,:)=points_three(i,:);
end
output=output';

%write the STL-file (without a loop)

STL_file = fopen (STL_file_name,'wt');

if (STL_file < 0)
    fprintf (1,'\n');
    fprintf (1,'Could not open the file "%s".\n',STL_file_name);
    error ('STL_WRITE - Fatal error!');
end

fprintf (STL_file,'solid %s\n',solid);
fprintf (STL_file, '  facet normal  %14e  %14e  %14e\n    outer loop\n      vertex %14e %14e %14e\n      vertex %14e %14e %14e\n      vertex %14e %14e %14e\n    endloop\n  endfacet\n',output);
   
time=toc;
fprintf (STL_file,'endsolid %s\n',solid );
fprintf (1,'\n');
fprintf (1,'attended time = %d\n',time);
fclose (STL_file);
%End of the function




