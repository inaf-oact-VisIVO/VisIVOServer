#version 120 
#extension GL_EXT_geometry_shader4 : enable


void main(void)
{
	int i;
    for(i=0; i< gl_VerticesIn; i++)
    {

        vec4 pos = gl_PositionIn[i]; 

        // Quad size will be 2*k
        // Multiply by 1000 to get around implementation defined point sprite size limit
        float k = gl_PointSizeIn[i] * 1000; 

        // Set colour
        vec4  c = vec4( gl_FrontColorIn[i].x, gl_FrontColorIn[i].y,gl_FrontColorIn[i].z, gl_FrontColorIn[i].w);

        gl_Position = pos + vec4(-k, -k,  0,  0);
        gl_TexCoord[0] =    vec4( 0,  0,  0,  0);
        gl_FrontColor = c;
        EmitVertex();

        gl_Position = pos + vec4( k, -k,  0,  0);
        gl_TexCoord[0] =    vec4( 1,  0,  0,  0);
        gl_FrontColor = c;
        EmitVertex();

        gl_Position = pos + vec4(-k,  k,  0,  0);
        gl_TexCoord[0] =    vec4( 0,  1,  0,  0);
        gl_FrontColor = c;
        EmitVertex();

        gl_Position = pos + vec4( k,  k,  0,  0);
        gl_TexCoord[0] =    vec4( 1,  1,  0,  0);
        gl_FrontColor = c;
        EmitVertex();

        EndPrimitive();
    }
}


