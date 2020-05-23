## Transverse field Ising mode: GS energies and Momentum matrices
I consider the Hamiltonian
$$
H = -\sum_{n=1}^{N-1}\sigma_n^x \sigma_{n+1}^x - \sigma_N^x \sigma_{1}^x - \sum_{n=1}^{N}\sigma_n^z.
$$
This is the Hamiltonian of the transverse field Ising model with periodic boundary conditions and transverse field set to $1$.
### Momentum matrix and GS energy for $N=2$
The ground state energy is:
$$
E_{GS}(N=2) = -2\sqrt{2}.
$$
Considering the vector

$$
\vec{v}^{\dagger} = (\mathbb{I}, \sigma^x_1,\sigma^x_2,\sigma^y_1,\sigma^y_2,\sigma^z_1,\sigma^z_2),
$$

the level $1$ momentum matrix $\Gamma := \langle v^{\dagger} v \rangle$  on the ground state is

$$
\Gamma = 
\begin{pmatrix}
1 	& 	0	&	0		&	0			&	0			&  \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}}\\
&1&\frac{1}{\sqrt{2}}&i\frac{1}{\sqrt{2}}&0&0&0\\
&&1&0&i\frac{1}{\sqrt{2}}&0&0\\
&&&1&-\frac{1}{\sqrt{2}}&0&0\\
&&&&1&0&0\\
&&&&&1&1&\\
&&&&&&1
\end{pmatrix}
$$

### Momentum matrix and GS energy for $N=3$
The ground state energy is:
$$
E_{GS}(N=3) = -4.
$$
Considering the vector

$$
\vec{v}^{\dagger} = (\mathbb{I}, \sigma^x_1,\sigma^x_2,\sigma^x_3,\sigma^y_1,\sigma^y_2,\sigma^y_3,\sigma^z_1,\sigma^z_2,\sigma^z_3),
$$

the level $1$ momentum matrix $\Gamma := \langle v^{\dagger} v \rangle$  on the ground state is

$$
\Gamma = 
\begin{pmatrix}
1 	& 	0	&	0		&	0			&	0			& 	0	&	0	&	2/3		&	2/3	&	2/3	 \\
 	& 	1	&2/3	&	2/3		&	i 2/3	& 	0	&	0	&	0			&		0		&	0			 \\
 	& 		&	1		&	2/3		&	0			& 	i 2/3	&	0			&		0		&	0			& 0 \\
	 & 		&			& 		1		&	0			&	0			& 	i 2/3	&	0			&		0		&	0	\\
&&&&1&-1/3&-1/3&0&0&0\\
&&&&&1&-1/3&0&0&0\\
&&&&&&1&0&0&0\\
&&&&&&&1&2/3&2/3\\
&&&&&&&&1&2/3\\
&&&&&&&&&1
\end{pmatrix}
$$

### Ground state energies
The ground state energy for a system of $N$ is:

$$
E_{GS}(N) = -\sum_{n=0}^{N-1} \sqrt{2+2\cos(\frac{\pi(1-N+2n)}{N})}
$$

|     $N$      |       $E_{GS}(N)$        |
| :--------: | :------------------------: |
|   $2$   |       $-2\sqrt{2}$      |
|   $3$    |           $-4$           |
|   $4$    | $-2\sqrt{2(2+\sqrt{2})}$ |
|   $5$    |     $-2(1+\sqrt{5})$     |
|   $6$    |  $-4\sqrt{2+\sqrt{3}}$   |
| $\vdots$ |         $\vdots$         |