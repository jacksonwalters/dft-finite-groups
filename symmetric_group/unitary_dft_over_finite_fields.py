{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "ba61a7c6-db85-4449-bd7a-b6209451b0a4",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute the uDFT by noting DFT.DFT^* = D, a diagonal matrix, and factoring as D = RR^*, so uDFT = R^{-1}.DFT\n",
    "#alternately, one could use unitary representations and normalization factors \\sqrt{d_\\rho/|G|}\n",
    "#the resulting DFT will have signs \\pm 1 on the diagonal DFT.DFT^* = S, which can be factored as s = rr^*, so S=RR^*\n",
    "#again, uDFT = R^{-1}.DFT"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "367e923c-c757-478e-bedd-e752e8792357",
   "metadata": {},
   "outputs": [],
   "source": [
    "#for u in GF(q), we can factor as u=aa^*=aa^q=a^{q+1} in GF(q**2) using gen. z and modular arithmetic\n",
    "def conj_sqrt(u):\n",
    "    if u == 0:\n",
    "        return 0\n",
    "    z = u.parent().multiplicative_generator()\n",
    "    k = u.log(z)  # Compute discrete log of u to the base z\n",
    "    if k % (q+1) != 0:\n",
    "        raise ValueError(\"exponent must be divisible by q+1\")\n",
    "    return z ** (k//(q+1))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "d436942b-2764-48cb-a318-e943960b671d",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute the uDFT by noting DFT.DFT^* = D, a diagonal matrix, and factoring as D = RR^*, so uDFT = R^{-1}.DFT\n",
    "def unitary_dft():\n",
    "    dft_matrix = SGA.dft()\n",
    "    sign_diag = (dft_matrix*dft_matrix.H).diagonal()\n",
    "    factor_diag_inv = diagonal_matrix([~conj_sqrt(d) for d in sign_diag])\n",
    "    return factor_diag_inv*dft_matrix"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "8638af45-82cd-4ee0-936d-3fdbeb482f32",
   "metadata": {},
   "outputs": [],
   "source": [
    "#parameters and define the symmetric group algebra\n",
    "n = 4; q = 7\n",
    "F = GF(q**2)\n",
    "SGA = SymmetricGroupAlgebra(F,n) # F[S_n], group algebra\n",
    "assert F.characteristic() > 0, \"F must have positive characteristic\"\n",
    "if not (F.is_field() and F.is_finite() and F.order().is_square()):\n",
    "    raise ValueError(\"the base ring must be a finite field of square order\")\n",
    "if F.characteristic().divides(SGA.group().cardinality()):\n",
    "    raise NotImplementedError(\"not implemented when p|n!; dimension of invariant forms may be greater than one\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "1251c71a-3b52-43e8-bf6d-8f3dea6e4007",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "24 x 24 dense matrix over Finite Field in z2 of size 7^2 (use the '.str()' method to see the entries)"
      ]
     },
     "execution_count": 5,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#converting the unitary DFT over finite fields to a complex matrix using the root of unity map\n",
    "U = unitary_dft(); U"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "58cdb5f2-a153-4884-9aab-906d8f9464b5",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "True"
      ]
     },
     "execution_count": 6,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#check that U is unitary over a finite field. .H is conjugate-transpose\n",
    "U*U.H == 1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "741f25cc-52c7-43c1-8d7c-7fc148df3bdc",
   "metadata": {},
   "outputs": [],
   "source": [
    "#given a multiplicative generator `z` of the finite field, the discrete_log is the exponent of the generator\n",
    "#the discrete_log of zero is -infinity, which we set to -1 for convenience since all other values are nonnegative\n",
    "discrete_log = lambda F, x: x.log(F.multiplicative_generator()) if x != 0 else -1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "0d5de1cd-1b8e-49b6-a959-0773c9c1314b",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute the discrete log of each entry of U\n",
    "log_U = U.apply_map(lambda x: discrete_log(F,x))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "a1624dfe-11e0-45fd-af51-7d5d70b86315",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "24 x 24 dense matrix over Integer Ring (use the '.str()' method to see the entries)"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "log_U"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "425df9fe-be30-40e3-8aff-fa1348bc74f4",
   "metadata": {},
   "outputs": [],
   "source": [
    "\"\"\"\n",
    "plot the discrete log of a matrix valued in a finite field\n",
    "since log(0) = -infinity, we set this to -1, and color it black\n",
    "\n",
    "- F is the field\n",
    "- M is a matrix of discrete log values of elements of F\n",
    "\n",
    "BUG: this is not working for large splitting fields, when the matrix values are ~7 digits\n",
    "\"\"\"\n",
    "def plot_discrete_log(F, M, path, title):\n",
    "    import numpy as np\n",
    "    import matplotlib.pyplot as plt\n",
    "    from matplotlib.colors import ListedColormap, BoundaryNorm\n",
    "\n",
    "    cmap = plt.cm.hsv  # Get the HSV colormap\n",
    "    num_colors = min(F.order(), 256) #for large fields we can't use that many colors, set a cutoff\n",
    "    new_colors = np.vstack(([0, 0, 0, 1], cmap(np.linspace(0, 1, num_colors))))\n",
    "    custom_cmap = ListedColormap(new_colors) # create a new custom colormap\n",
    "    norm = BoundaryNorm([-1] + list(np.linspace(0, F.order()-1, num_colors)), custom_cmap.N) #map -1 to black\n",
    "\n",
    "    # Plotting the data\n",
    "    plt.imshow(M, cmap=custom_cmap, norm=norm, interpolation=\"nearest\")\n",
    "    plt.title(title, fontsize=16)\n",
    "    plt.colorbar()\n",
    "    plot_title = path + title.replace(' ','_') + '.png'\n",
    "    plt.savefig(plot_title, dpi=300, bbox_inches=\"tight\")\n",
    "    plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "dd212ab8-22f9-4968-954f-64a822d3e205",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAecAAAG3CAYAAACOg52jAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjguMCwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy81sbWrAAAACXBIWXMAAA9hAAAPYQGoP6dpAABD9klEQVR4nO3df3gTVb4/8PdQ2vQHTbDUtimUWLTIjwJKi6UFoaCtdr1coKgoLlv2WREvP7T26+IK6xJdaF3uXcRrFRd0Ea6wsHsF9D4oUJUWkWW3IKxdcBW0QFkplQIJ/ZVCe75/YGYNSdskJyWT9v16njxtZubM58zJJJ9MZs4cRQghQERERJrRw98VICIiIkdMzkRERBrD5ExERKQxTM5EREQaw+RMRESkMUzOREREGsPkTEREpDFMzkRERBrD5ExERKQxAZ+cb7rpJiiKghMnTjhMz8zMhKIoKC0t9Uu9uoqu3o6lpaWYMGEC9Ho9FEVxuS+Ra7W1tZg9ezb69u2LoKAgKIoCs9ns72oRdQkBn5y7ixMnTsBsNuOtt97yd1W6jCNHjuCee+5BaWkpoqOjMWbMGIwZMwahoaGdFtP+BcD+CAoKQu/evXHTTTfh3/7t37Bs2TKcOnWq3XXYv5C291i5ciVmzZrV4XKuHu5+OZk8eTLeeOMN1NfXIzU1FWPGjEH//v190EpdixAC48aNU9t37969/q6Sx9566y239x/yjZ7+rkBn6d+/P2699VaEh4f7uyo+ceLECTz//PMYP348Zs2a5e/qdAlvvvkmmpubsWDBAvz3f//3dY2dnJwMg8EAAGhoaEBNTQ22b9+O7du341e/+hUee+wx/Pa3v213/01KSkJMTIzLeX379kVDQwPGjBnjNO/AgQOw2Wxtlnfny8nnn3+OTz/9FH379sWRI0fUbSFnb775Jj755BN/V0NKbGysy33J7h//+Adqa2uRnp5+HWvVtXXZ5Lx+/Xp/V4E07h//+AcAICcn57rHfuWVV5CZmekw7Z///Cd+//vf48UXX8Trr7+Oo0eP4sMPP0RwcLDLdSxatKjDL2qLFi1ymnbTTTfh5MmTbpVvi73txowZw8Tcju+++w7PPPMMbr/9dnz33Xc4ffq0v6vklZycnDbfJ0IImEwm1NbWYubMmde5Zl0Xf9ambquxsREAEBYW5ueaXNW3b18899xz2Lt3LyIiIrBnzx4UFhb6u1ouaa3ttOqpp57ChQsX8NprryEoKMjf1ekUpaWlqKqqQkhICKZPn+7v6nQZAZGcT548iR//+MeIiYlBeHg4hg8fjldffRXtjXbZ1oVMV65cwcsvv4w77rgDkZGR0Ol0iI+PR0ZGBpYsWYKLFy86revKlStYs2YNJkyYgD59+iA0NBQDBgzAtGnT8O6777YZ9/Dhw7j//vsRGxuLHj16OJwvvnLlCl5//XWMHTsWvXv3RmhoKAYNGoRf/vKXsFqtTuucMGECAKCsrMzh/M5NN93kVN+//vWveOihh9C3b1+EhIQgNjYWDzzwAA4dOtR+Q3tICIG3334b48ePR+/evREWFoZBgwbhmWeewfnz59ss9/XXX+Phhx/GjTfeiPDwcNx22214/fXXAbR9gV9HLl++jFdeeQV33HEH9Ho9IiIiMGLECCxbtgwNDQ0Oy9rPxdr3jQkTJqjt6c6RZGlpKRRFcTrytTtx4kSbr407br/9dixduhQA8PLLLzvV35/s225vp3Xr1rV5vrG+vh5Lly7F8OHDERERAb1ej7S0NLz66qu4cuVKm+vOzMzElStXsHz5cgwbNgzh4eFut+UP95/9+/cjJycHN9xwAyIiInDnnXfi448/lm0Ct3344YfYsGEDHn30UYwePVp6fe19Dnr7vvGF//mf/wEA/OhHP0JUVNR1j99lCY07evSo6NOnjwAgQkNDRUpKiujfv78AIObOnStMJpMAICorKx3KjR8/XgAQu3fvdpg+bdo0AUAAEDfffLMYNWqUSEhIEEFBQQKAOHTokMPy58+fF2PGjFHLmEwmkZqaKmJiYtTnruI+//zzQqfTiV69eomUlBQxYMAAsXbtWiGEEBaLRYwbN04AED169BAmk0kkJyeLkJAQAUAMHjxYnD17Vl3n/PnzRXJysgAg9Hq9GDNmjPq4//77HeKvWLFCKIoiAIioqChx++23q+0XHBws3nnnHY/av612bG1tFTNmzFDbZcCAAWLkyJHqNphMJvH11187re9vf/ub6N27twAgwsLCREpKivoaPvHEE22+nu1paGgQEydOVOsyePBgMXz4cNGjRw8BQNx2223i3Llz6vLLli0TY8aMEXq9XgAQycnJansuW7asw3i7d+8WAMT48eNdzq+srHS5bwgh1Dpe257XslqtIjg4WAAQH374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",
      "text/plain": [
       "<Figure size 640x480 with 2 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plot_discrete_log(F,log_U,path='plots/dft_matrix/',title=f\"discrete log of uDFT for n={n} q={q}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "dbf3c83d-ea6e-4dcd-a491-332a715dcbe1",
   "metadata": {},
   "outputs": [],
   "source": [
    "#the map from modular representation theory to compute Brauer characters from F_q --> \\C giving roots of unity\n",
    "#note: Brauer character is the rep'n matrix eigenvalues (over a splitting field extension of F_q) mapped to \\C using Brauer map and summed\n",
    "#i.e. let \\alpha = g^k |--> \\exp(2*pi*i*k/(q-1))\n",
    "#i.e. F_q^* is cyclic of order q-1, mapping to (q-1)^th roots of unity in \\C\n",
    "#if the discrete log is provided, we use it directly; otherwise, we compute it\n",
    "def brauer_map(F, a=None, log_a=None):\n",
    "    \"\"\"\n",
    "    Map from F_q to C using the Brauer character formula.\n",
    "    If `log_a` is provided, it uses that directly; otherwise, it computes the log.\n",
    "    \n",
    "    a: Element of F_q or None if log_a is given\n",
    "    log_a: Precomputed log value (optional)\n",
    "    F: Finite field F_q\n",
    "    \n",
    "    Returns: Complex value of the Brauer character\n",
    "    \"\"\"\n",
    "    if a is None and log_a is None:\n",
    "        raise ValueError(\"Either 'a' or 'log_a' must be provided.\")\n",
    "    \n",
    "    if a is not None:\n",
    "        if a == 0:\n",
    "            return 0\n",
    "        log_a = a.log(F.multiplicative_generator())\n",
    "    \n",
    "    return exp(2 * pi * I * log_a / (F.order() - 1))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "913954c0-dfdf-4099-b9b2-7c59a5e6e767",
   "metadata": {},
   "outputs": [],
   "source": [
    "#complexify the uDFT matrix over F_q using the Brauer map\n",
    "U_complex = matrix(CC,U.apply_map(lambda a: brauer_map(F,a)))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "2ce4d9b0-c891-40de-9e2d-92e941e05901",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute the Gram matrix, taking inner products of rows and columns w.r.t conjugate inner product\n",
    "#note: want this to be the identity matrix (so it would be unitary over \\C), but currently it is not quite\n",
    "#unitary matrices are required if they are to be used as operators in quantum computing\n",
    "gram = U_complex*U_complex.H"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "9e4e0a10-8c5f-42df-974b-0e729c6ba084",
   "metadata": {},
   "outputs": [],
   "source": [
    "#function to round each component of a complex number \n",
    "def round_complex(z, digits):\n",
    "    if z.imag_part():\n",
    "        return round(z.real_part(), digits) + round(z.imag_part(), digits) * I\n",
    "    return round(z, digits)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "2ce34d19-f35e-486b-b373-1c77c4e5724f",
   "metadata": {},
   "outputs": [],
   "source": [
    "#round the (complex) Gram matrix to three decimal places\n",
    "gram_rounded = gram.apply_map(lambda u:round_complex(u,3))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "cc3df8b2-4ee4-4643-8963-be7578c9485c",
   "metadata": {},
   "outputs": [],
   "source": [
    "#plot the complexified version of the uDFT matrix over a finite field\n",
    "def plot_arg_complex(U_complex, title):\n",
    "    U_arg = U_complex.apply_map(lambda x: arg(x))  # find the argument of each element\n",
    "    plot = matrix_plot(U_arg, cmap='hsv', colorbar=True, title=title)  # plot the matrix\n",
    "    filename = \"plots/dft_matrix/\" + title.replace(\" \", \"_\") + \".png\"\n",
    "    plot.save(filename, dpi=300)  # Save the plot as a PNG file with high resolution\n",
    "    return plot"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "daca587f-859b-42fb-9c22-08f9e377c953",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAhIAAAHWCAYAAAAirGCAAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjguMCwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy81sbWrAAAACXBIWXMAAA9hAAAPYQGoP6dpAAA/AklEQVR4nO3de3xU1bn/8e8QSAiQBCOSi4SAoEiMggJV8BJoFclpEcUqCirUy0/KYKVobRFbIlajFjnaMmLVFtCCWOu1ahG8ELwcNSBUahCBRolKmiPaBFHJbf3+oMxhTIBk1prMnszn/XrtF8xlrfXM3nsyzzxrz94+Y4wRAABAGDpEOwAAABC7SCQAAEDYSCQAAEDYSCQAAEDYSCQAAEDYSCQAAEDYSCQAAEDYSCQAAEDYSCQAAEDYSCSgr776Sueff75SU1Pl8/n073//O9ohtYrP59NTTz3lrL8pU6bo3HPPDd42xuj//b//p/T0dPl8Pm3YsEEjR47UjBkzrMZZvHixunfvbtVHLHjqqafUv39/JSQkWK8zAN5DIgEtWbJEr776qt544w3t2LFDaWlp0Q4pqu655x4tXrw4eHvFihVavHixnn32We3YsUP5+fl64okndMstt0QvyP9YvXq1fD6ffD6fOnTooLS0NJ144om64YYbtGPHjpDnFhUVBZ+7//Lggw82e//+S1FRUdgxXn311frhD3+oioqKA66z9evX6wc/+IF69uypzp07q0+fPpowYYI+++yzVo21c+dO9erVK+YS4gNtmxdffLFV/Sxfvlw+ny8kEZakXbt2acaMGcrNzVVycrJGjBih0tJSh68A8axjtANAy9TV1alTp04R6Xvbtm0aOHCg8vPzI9J/rPl2IrVt2zZlZWVpxIgRwfvS09PbOqyD2rx5s1JTU1VTU6N33nlHd955p/7whz9o9erVOv7444PPO+6445p8OB122GH6wQ9+ELw9b948rVixIuR53bp1CyuuL7/8UlVVVTr77LOVnZ3d7HOqqqp05plnauzYsXrhhRfUvXt3lZeX65lnntFXX33VqvGuuOIKnXDCCfrkk0/CijfSamtrlZiY2OxjzW2b1uxnH330ka6//nqdfvrpTR678sor9Y9//EMPP/ywsrOz9ac//UlnnnmmysrKdOSRR7buRQDfZtDm/va3v5lTTz3VpKWlmfT0dPP973/fbN26Nfh4eXm5kWQeffRRU1BQYJKSkswf//hHU1dXZ6655ppguxtuuMFcdtllZty4cQcd7y9/+YvJy8sziYmJJjc318ybNy/4WEFBgZEUXAoKCg7Yz9NPP22GDBlikpKSzOGHH27OO++84GOff/65ufTSS0337t1NcnKyGTNmjPnggw+Cjy9atMikpaWZv/71r+aYY44xycnJ5vzzzzdffvmlWbx4scnNzTXdu3c306dPN/X19cF2ubm5Zu7cuebiiy82Xbt2NVlZWea3v/1tSFySzJNPPhm8/fHHH5sLL7zQdO/e3aSnp5tzzjnHlJeXG2OM2bRpk0lOTjZLly4NPv/xxx83SUlJ5t133zXGGDN58uTgOp08eXLI+snNzQ2ut2uvvTbYx549e8zPfvYzk52dbbp06WK+853vmFdeeSUkzkWLFpmcnByTnJxszj33XDNv3jyTlpZ2wPX9yiuvGEnmiy++CN63fv16Iyn4epp7jjHGfPXVV2bAgAHm1FNPDd43Z84cM2jQoAOO19rnGXPw7b4vtv2Xb68TY4x58sknTceOHU1dXV2LxjyQe++91xQUFJiXXnqp2XXybR999JE555xzTNeuXU1KSoq54IILTGVlpTHGmPfff99IMps2bQppc9ddd5nc3FzT2NhojDHmvffeM4WFhaZr166mZ8+e5pJLLjH/+7//G3x+QUGB8fv95qc//ak5/PDDzRlnnNFsLK1Z582pr683p556qnnwwQdD9l9j9u4LCQkJ5tlnnw1pM2jQIDN79uywxwT2YWojCnbv3q2ZM2eqtLRUL730kjp06KDzzjtPjY2NIc/7+c9/rp/85CfatGmTzj77bN1xxx1aunSpFi1apNdff101NTWHPDZg3bp1uvDCC3XRRRdp48aNKioq0i9/+ctg6f6JJ57QVVddpeHDh2vHjh164oknmu3nueee0/jx4/X9739f69ev10svvaShQ4cGH58yZYrWrl2rZ555Rv/zP/8jY4z+67/+S3V1dcHnfPXVV/rtb3+r5cuXa8WKFVq9erXGjx+v559/Xs8//7wefvhh3X///frLX/4SMvZvfvMbnXDCCXrnnXc0a9Ys/fSnP9WqVauajfOrr77SqFGj1K1bN61Zs0avvfaaunXrpjFjxqi2tlbHHnus5s2bp2nTpumjjz7Sp59+qquuukq33357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YyeddJIef/xx9enTRx07RudPaXNVm5tuukm7du3SPffcE1JxkfYmNV27dtUXX3yhF154QXfeeWdbhot2iqmNNnbYYYfp8MMP1/3336+tW7fq5Zdf1syZM1vU9pprrlFxcbGefvppbd68Wddee62++OKLg377u+666/TSSy/plltu0QcffKAlS5ZowYIFuv7661sV95w5c/TII49ozpw52rRpkzZu3Bj8I3T00Udr3Lhxuuqqq/Taa6/p73//uy655BIdeeSRGjduXKvGac7rr7+uO++8Ux988IECgYAee+wxXXvttc0+d9KkSerRo4fGjRunV199VeXl5SopKdG1116rjz/+WJI0depU5eTk6KabbtL8+fNljGn1+tjfMccco0mTJumyyy7TE088ofLycpWWluqOO+7Q888/L0n6yU9+ohUrVgRfx4IFC7RixYqD9tu/f3/l5OSoqKhIH3zwgZ577jndddddzT63qqpKlZWV2rJli5YvX65TTz1Vn332mRYuXBj26zoUV9v92Wef1SWXXKJnn31WH3zwgTZv3qx58+bp+eefb1E//fr1U35+fnDZ9y174MCB6tmzZ7NtzjzzTJ1wwgnBJPLtt9/WZZddpoKCgpCpiPHjx6umpkY//vGPNWrUqJCfSvr9fn3++ee6+OKL9fbbb+uf//ynVq5cqcsvv7xV04Y2OnfuHPLa8/Pz1b17d6WkpCg/Pz/4U9MXXnhBK1asUHl5uVatWqVRo0ZpwIABzX4xAFqLRKKNdejQQcuXL9e6deuUn5+vn/70p/rNb37Torb7vhlddtllGj58uLp166azzz67ydHY+zvppJP05z//WcuXL1d+fr5+9atfae7cuZoyZUqr4h45cqQee+wxPfPMMxo8eLC++93v6q233go+vmjRIg0ZMkQ/+MEPNHz4cBlj9Pzzzzs598V1112ndevW6cQTT9Qtt9yiu+66S2effXazz+3SpYvWrFmj3r17a/z48Ro4cKAuv/xyff3110pNTdVDDz0U/IVIx44d1aVLFy1dulQPPvhg8EM/HIsWLdJll12m6667TgMGDNA555yjt956K/iN8JRTTtGDDz6o3/3udxo8eLBWrlypm2666aB9durUSY888ojef/99DRo0SHfccccBf+kxYMAAZWdna8iQIbr99tt15pln6h//+EdwPj9SXGz3vLw8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",
      "text/plain": [
       "Graphics object consisting of 1 graphics primitive"
      ]
     },
     "execution_count": 18,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#plot the complexified version of the uDFT matrix over a finite field\n",
    "plot_arg_complex(U_complex,title=f\"arg of complexified uDFT of S_{n} over F_{q**2}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "1e1f300a-33e1-435d-9e15-8548b27ea103",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "True"
      ]
     },
     "execution_count": 19,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#check that Gram matrix is conjugate symmetric (should be since it is U*U.H)\n",
    "gram_rounded == gram_rounded.H"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "id": "89edaab1-11b3-4425-b764-25a301e78785",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAhIAAAHWCAYAAAAirGCAAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjguMCwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy81sbWrAAAACXBIWXMAAA9hAAAPYQGoP6dpAABCXElEQVR4nO3deXxU1f3/8fckJgFCEglLFg0hKlQCSBUigiJQJZK2VOq+FENdKhoQpCpulbgRd211QLEKWBesittXq+ACqOiXoOSnEoSgKKkQqaAJBM16fn/wzZQhCSRzz2TuJK/n43EfOnfuOfczd26Yz3zOuXc8xhgjAACAAESEOgAAABC+SCQAAEDASCQAAEDASCQAAEDASCQAAEDASCQAAEDASCQAAEDASCQAAEDASCQAAEDASCRaaPfu3Tr99NMVHx8vj8ejH3/8MdQhhTWPx6OXXnopqPsoKyvT2LFjFRsbq4MPPjio+7Lt66+/lsfjUVFRkbU+R48erenTp/seN3VO9+nTRw888ICj/eTn5+uXv/yloz7Cwbx585SWlqaIiAjHxwwIZyQSLbRw4UK99957WrlypbZu3aqEhIRQh+Qqy5Yta1WCtXXrVuXk5AQ1pvvvv19bt25VUVGRNmzYENR9hYPFixfr1ltv9T1u6pwuLCzUn/70pxBGuceCBQvk8Xjk8XgUGRmpbt26adiwYbrllltUXl7ut+2kSZN82+69/P3vf29y/d7LggULAoqvoqJCU6ZM0cyZM/Xtt982e8zeffddjRkzRomJierSpYv69u2r3Nxc1dbWtmp/GzduVFxcXNglxM29Nxs3bmxVPwUFBfJ4PH6JsCR99913mjRpklJTU9WlSxeNGzdOJSUlFl8BWuKgUAdgU01NjaKiooLS95dffqn+/ftr4MCBVvutrq5WdHS01T7drOH1JicnB31fX375pYYMGaK+ffsGfV/hIDEx0e9xU+d0z5492zqsZsXHx2v9+vUyxujHH3/UypUrVVBQoPnz5+uDDz5Qamqqb9tx48Zp/vz5fu27deum3/72t77H06ZNU0VFhd92gX4h2Lx5s2pqavSb3/xGKSkpTW6zdu1a5eTk6IorrtCDDz6ozp07q6SkRM8//7zq6+tbvK+amhqde+65GjlypFauXBlQvMG2v3/HmnpvWnOeFRYWat68eTrqqKP81htjNGHCBEVFRenll19WfHy87rvvPp188skqLi5WbGxs618IAmNc6l//+pc5/vjjTUJCgklMTDS/+c1vzMaNG33Pb9q0yUgyzz77rBk1apSJiYkxjz/+uKmpqTFTp071tbvmmmvMBRdcYE499dT97u/55583mZmZJjo62qSnp5t77rnH99yoUaOMJN8yatSoZvu59dZbTc+ePU3Xrl3NRRddZGbOnGkGDx7sez43N9eceuqpZvbs2SYlJcWkp6cbY4z5xz/+YYYMGWK6du1qkpKSzLnnnmu+++47X7t3333XSDJvvPGG+eUvf2k6depkxowZY7777jvz+uuvmyOPPNLExcWZc845x1RWVjYb3/z5801CQoJ59dVXTb9+/Uznzp3N6aefbnbt2mUWLFhg0tPTzcEHH2ymTJliamtrfe32F1/De7H3kpub6zt2eXl55sorrzTdu3c3J554ojHGGEnmxRdfNMYYs3DhQhMbG2s2bNjg29+UKVNM3759za5du5p9LXPmzDGHHXaYiYqKMv369TNPPPGE77n09PQm42nKY4895nvvk5OTTV5enu+5b775xvzud78zsbGxJi4uzpx55pmmrKzM9/ysWbPM4MGDzWOPPWbS0tJMbGysmTx5sqmtrTV33nmnSUpKMj179jS33Xab3z4lmTlz5phx48aZTp06mT59+ph//vOfvucbjumaNWt869auXWtycnJMbGys6dWrl/nDH/5g/vOf/xhj9pwfUVFRZsWKFb7t77nnHtO9e3ezZcsW33sxbdo03/83dU6np6eb+++/39fHjz/+aC655BLTs2dPExcXZ8aMGWOKior8XktBQYHp1auX6dq1q7nwwgsbnfP7ajgH9/biiy+avf85amobY4z57rvvTI8ePcz555/vW9fwN3UgLd3OmP2/7/Pnz290vm/atKlRH/fff7/p06dPi/a3P9dcc435wx/+0Owx2denn35qxowZYzp16mQSExPNJZdcYnbu3GmMMeaNN94wMTEx5ocffvBrM3XqVN/fpjHGfPDBB2bkyJGmU6dO5tBDDzVTp071+1tMT083t956q8nNzTXx8fHmggsuaDKW1hzzpuzcudP07dvXLF261O/8NcaY9evXG0nm888/962rra01iYmJ5tFHHw14n2g91yYSzz//vHnhhRfMhg0bzJo1a8z48ePNoEGDTF1dnTHmv//Q9unTx7zwwgvmq6++Mt9++6257bbbTGJiolm8eLFZt26dmTx5somPj9/vybx69WoTERFhbrnlFrN+/Xozf/5807lzZzN//nxjjDHbt283l1xyiRk+fLjZunWr2b59e5P9PPnkk6ZTp07m8ccfN+vXrzc333yziY+Pb5RIdO3a1UycONF8/vnn5rPPPjPG7Pkwe/31182XX35pPvzwQ3PccceZnJwcX7uGROK4444z77//vvnkk0/MEUccYUaNGmWys7PNJ598YlasWGG6d+9u7rjjjmZf6/z5801UVJQZO3as+eSTT8zy5ctN9+7dTXZ2tjnrrLPM2rVrzauvvmqio6PNokWLfO32F19tba154YUXjCSzfv16s3XrVvPjjz8aY/Z8YHXt2tVcffXV5osvvjDr1q0zxvgnEsYYc+aZZ5qsrCxTU1Nj/vWvf5moqCizatWqZl/H4sWLTVRUlPF6vWb9+vXm3nvvNZGRkeadd94xxhizbds2M27cOHPWWWf5xbOvOXPmmE6dOpkHHnjArF+/3qxatcr3QVpfX2+OPvpoc8IJJ5jVq1ebjz76yBxzzDF+ieSsWbNM165dzRlnnGHWrl1rXnnlFRMdHW1OOeUUM3XqVPPFF1+Yxx9/3EgyH374oa+dJNO9e3fz6KOPmvXr15sbb7zRREZGmuLiYmNM40Riy5YtpkePHua6664z69atM5988okZO3asGTNmjK/Pq6++2qSnp5sff/zRFBUVmZiYGLN48WLf83v/Q9zcOb13IlFfX2+OP/54M378eFNYWGg2bNhg/vznP5vu3bv7tn/22WdNdHS0efTRR80XX3xhbrjhBhMXFxe0RMIYY6ZNm2bi4uJ8ia7tROJA7/vu3bvNW2+9ZSSZVatWma1bt/ol3Q2eeeYZExMTY5YvX37AfTbn7bffNhkZGaa8vLxFiURlZaVJTU01p512mvnss8987RsS6draWpOUlGT+/ve/+9o0rHvkkUeMMXsSka5du5r777/fbNiwwXzwwQfm6KOPNpMmTfK1SU9PN/Hx8ebuu+82JSUlpqSkpMl4nCYSF1xwgZk+fboxxjRKJD799FMjye8LpjHGJCcn7/eLA+xzbSKxr23bthlJvg/ehn9oH3jgAb/tkpKSzN133+17XFtba3r37r3fk/m8884zY8eO9Vt39dVXm8zMTN/jadOm7bcSYYwxw4YN8/s2a4wxxx9/fKNEIikpyVRVVe23r1WrVhlJvm8SDYnEW2+95dumoKDASDJffvmlb92ll15qTjnllGb7bfg2tfcf36WXXmq6dOni25cxxpxyyinm0ksvbXV8+37TGTVqlPnlL3/ZqP2+icSOHTvMoYceai677DKTlJTU6Bv8vkaMGGEuueQSv3Vnnnmm+fWvf+17fOqppx7wH5TU1FRzww03NPnckiVLTGRkpNm8ebNv3dq1a30fIMbsSSS6dOliKioqfNuccsoppk+fPr6k1xhjfvGLX5iCggLfY0lm8uTJfvsbNmyYueyyy4wxjROJv/zlLyY7O9tv+9LSUl/yZowxVVVV5uijjzZnnXWWGTBggLn44ov9tt/3H+Kmzum9E4m3337bxMfHm59//tlvm8MPP9z3oTN8+PAmX0cwE4m5c+caSb6KWG5uromMjDSxsbG+5YwzzmjUrqUfai1539esWdNsJaJBbW2tmTRpkpFkkpOTzYQJE8yDDz5oysvLDxiDMcZ8//33Ji0tzZeItCSRmDdvnunWrZtf9eC1114zERERvorKFVdcYX71q1/5nn/zzTdNdHS02bFjhzHGmIkTJ5o//elPfv2+9957JiIiwvz000/GmD3nyYQJEw74Glr63jTlmWeeMQMHDvTtc9/zt7q62qSnp5szzzzT7Nixw1RVVfn+Tdz3bwXB5drJll9++aXOO+88HXbYYYqPj1dGRoakPWOTexs6dKjv/8vLy/Xdd9/p2GOP9a2LjIzUkCFD9ruvdevW6fjjj/dbd/zxx6ukpER1dXUtjnn9+vV++5bU6LEkDRo0qNF44po1a3TqqacqPT1dcXFxGj16tKTGr3fvccKkpCR16dJFhx12mN+6bdu27TfOLl266PDDD/dr06dPH3Xt2rXZfloaX1P2fo+a061bNz322GOaO3euDj/8cF177bX73b6592zdunUH3FeDbdu2acuWLTrppJOa3UdaWprS0tJ86zIzM3XwwQf77adPnz6Ki4vzPU5KSlJmZqYiIiL81u37vgwfPrzR4+bi//jjj/Xuu++qa9euvuXII4+UtOdvRZKio6P15JNP6oUXXtBPP/3k+EqCjz/+WLt27VL37t399rtp0ybfPtetW9fk6wgmY4ykPVf+NBgzZoyKiop8y9/+9reA+2/p+34gkZGRmj9/vv7973/rrrvuUmpqqm6//XYNGDBAW7duPWD7Sy65ROedd55OPPHEVsU+ePBgv/kBxx9/vOrr67V+/XpJ0vnnn69ly5Zpy5YtkqSnnnpKv/71r9WtWzdJe973BQsW+L3np5xyiurr67Vp0yZfvy35u5YCe29KS0s1bdo0Pfnkk+rUqVOT20RFRemFF17Qhg0bfJNZly1bppycHEVGRrYoNtjh2smW48ePV1pamh599FGlpqaqvr5eAwcOVHV1td92TU2o2fsfGOm///A0xxjT6jbNaUk/+8ZcWVmp7OxsZWdn68knn1TPnj21efNmnXLKKY1e796TST0eT6PJpR6P54ATuZpqs79+WhNfU1o66WnFihWKjIzUli1bVFlZqfj4+P1u39Sx3nfd/nTu3Hm/zzfX377rW3s896e5+Ovr6zV+/HjdeeedjZ7be7Jfw2S8HTt2aMeOHY4mnNXX1yslJUXLli1r9JyTqwciIiIa/V3U1NS0uP26desUHx+v7t27+9bFxsbqiCOOCDimvbX0fW+pQw45RBMnTtTEiRN12223qV+/fnr44Yd1880377fdO++8o1deeUX33HOPb//19fU66KCDNG/ePF144YWtirFh/bHHHqvDDz9cixYt0mWXXaYXX3zRbzJkfX29Lr30Ul1xxRWN+ujdu7fv/1t6bgXy3nz88cfatm2b35fAuro6rVixQg899JCqqqp8XxKLiopUXl6u6upq9ezZU8OGDWtxkgM7XFmR2L59u9atW6cbb7xRJ510kvr3768ffvjhgO0SEhKUlJSkVatW+dbV1dVpzZo1+22XmZmp999/32/dypUr1a9fv1Zltr/4xS/89i1Jq1evPmC7L774Qt9//73uuOMOjRw5UkceeeQBqwptqSXxNVRYWlPB2dvKlSt111136dVXX1V8fLymTp263+379+/f5HvWv3//Fu8zLi5Offr00dtvv93k85mZmdq8ebNKS0t964qLi1VeXt6q/TTno48+avS4ocqwr2OOOUZr165Vnz59dMQRR/gtDf+gf/nll7ryyiv16KOP6rjjjtMFF1zQqqsDmtpnWVmZDjrooEb77NGjh6Q970NTr2N/evbsqZ07d6qystK3rqX3y9i2bZuefvppTZgwwa/iY1Mw3/du3bopJSXF77U358MPP/T7Jn/LLbcoLi5ORUVF+v3vf99s7EVFRX79f/DBB4qIiFC/fv1868477zw99dRTevXVVxUREaHf/OY3vucazrV93/Mjjjiiza4wO+mkk/TZZ5/5vf6hQ4fq/PPPV1FRUaN/lxMSEtSzZ0+VlJRo9erVOvXUU9skTuzhykSiW7du6t69u+bNm6eNGzfqnXfe0YwZM1rUdurUqSooKNDLL7+s9evXa9q0afrhhx/2+03iz3/+s95++23deuut2rBhgxYuXKiHHnpIV111Vavinjp1qh577DEtXLhQJSUluu222/Tpp58e8FtM7969FR0drQcffFBfffWVXnnlFb/r/UOtJfGlp6fL4/Hof/7nf/Sf//xHu3btanH/O3fu1MSJEzV16lTl5OTo6aef1j//+U8999xzzba5+uqrtWDBAj388MMqKSnRfffdp8WLF7f6PcvPz9e9996rv/3tbyopKdEnn3yiBx98UJJ08skn66ijjtL555+vTz75RKtWrdIFF1ygUaNGWfnG89xzz+nxxx/Xhg0bNGvWLK1atUpTpkxpctu8vDzt2LFD5557rlatWqWvvvpKS5Ys0YUXXqi6ujrV1dVp4sSJys7O1h//+EfNnz9fn3/+ue69996A4zv55JM1fPhwTZgwQW+++aa+/vprrVy5UjfeeKMvQZ42bZoef/xxv9exdu3a/fY7bNgwdenSRddff702btyop59+usn7ORhjVFZWpq1bt2rdunV6/PHHNWLECCUkJOiOO+4I+HUdiK33/ZFHHtFll12mJUuW6Msvv9TatWs1c+ZMrV27VuPHjz9g+4ZLcxuWQw45RBERERo4cKBvGGJf559/vjp16qTc3Fx9/vnnevfddzV16lRNnDhRSUlJftt98sknuv3223XGGWf4DR/MnDlTH374ofLy8lRUVKSSkhK98sorB0zubYqLi/N77QMHDlRsbKy6d+/ud7nyc889p2XLlumrr77Syy+/rLFjx2rChAnKzs5us1jh0kQiIiJCixYt0scff6yBAwfqyiuv1N13392itjNnztS5556rCy64QMOHD/eN7zU3zibtycD/+c9/atGiRRo4cKBuuukm3XLLLZo0aVKr4j7//PN13XXX6aqrrtIxxxyjTZs2adKkSfvdt7TnG9qCBQv03HPPKTMzU3fccYevnOkGLYnvkEMO0c0336xrr71WSUlJzX4gNmXatGmKjY3V7NmzJUkDBgzQnXfeqcmTJ+vbb79tss2ECRP017/+VXfffbcGDBigRx55RPPnz/fN3Wip3NxcPfDAA5ozZ44GDBig3/72t74b2jTcfbNbt2468cQTdfLJJ+uwww7Ts88+26p9NOfmm2/WokWLdNRRR2nhwoV66qmnlJmZ2eS2qamp+uCDD1RXV6dTTjlFAwcO1LRp05SQkKCIiAjdfvvt+vrrrzVv3jxJUnJysv7+97/rxhtvDPjumB6PR6+//rpOPPFEXXjhherXr5/OOeccff31174PpbPPPls33XSTZs6cqSFDhuibb77RZZddtt9+ExMT9eSTT+r111/XoEGD9Mwzzyg/P7/RdhUVFUpJSdEhhxyi4cOH65FHHlFubq7WrFnT7L0bbLD1vh977LHatWuXJk+erAEDBmjUqFH66KOP9NJLL2nUqFFBib1Lly568803tWPHDmVlZemMM87QSSedpIceeshvu759+yorK0uffvqpzj//fL/njjrqKC1fvlwlJSUaOXKkjj76aP3lL38J6jEP1NatWzVx4kQdeeSRuuKKKzRx4kQ988wzoQ6rw/GYQCcDhIn6+nr1799fZ511Vki+5Y8dO1bJycn6xz/+0eb7hnt5PB69+OKLmjBhQqhDAQBHXDvZMlDffPONlixZolGjRqmqqkoPPfSQNm3apPPOOy/o+969e7cefvhhnXLKKYqMjNQzzzyjt956S0uXLg36vgEACAVXDm04ERERoQULFigrK0vHH3+8PvvsM7311ltWJscdSEMpeOTIkRoyZIheffVVvfDCCzr55JODvm8A4SUnJ8fvEsu9l4ZhvvZq8+bNzb72rl27tuiycrhHux/aAAA3+vbbb/XTTz81+VxiYmKj30ZpT2pra/X11183+3yfPn100EHtrmDebpFIAACAgLW7oQ0AANB2SCQAAEDASCQAAEDASCQAAEDA2l0iMWfOHGVkZKhTp04aMmSI3nvvvVCHFBby8/Pl8Xj8luTk5FCH5WorVqzQ+PHjlZqa6rsb4t6MMcrPz1dqaqo6d+6s0aNHH/D20R3NgY7hpEmTGp2Xxx13XGiCdaGCggJlZWUpLi5OvXr10oQJE3y/8tmA8xDB1q4SiWeffVbTp0/XDTfcoDVr1mjkyJHKycnhmuQWavh544bls88+C3VIrlZZWanBgwc3uv1wg7vuukv33XefHnroIRUWFio5OVljx47Vzp072zhS9zrQMZSkcePG+Z2Xr7/+ehtG6G7Lly9XXl6ePvroIy1dulS1tbXKzs72+9EuzkMEnWlHjj32WDN58mS/dUceeaS59tprQxRR+Jg1a5YZPHhwqMMIW5LMiy++6HtcX19vkpOTzR133OFb9/PPP5uEhATz8MMPhyBC99v3GBpjTG5urjn11FNDEk842rZtm5Fkli9fbozhPETbaDcVierqan388ceNfvUtOztbK1euDFFU4aWkpESpqanKyMjQOeeco6+++irUIYWtTZs2qayszO98jImJ0ahRozgfW2nZsmXq1auX+vXrp0suuaTRT9jjv8rLyyXJdzMrzkO0hXaTSHz//feqq6vz+6lcSUpKSlJZWVmIogofw4YN0xNPPKE333xTjz76qMrKyjRixAht37491KGFpYZzjvPRmZycHD311FN65513dO+996qwsFC/+tWvVFVVFerQXMcYoxkzZuiEE07w/dQ25yHaQru7B6nH4/F7bIxptA6N5eTk+P5/0KBBGj58uA4//HAtXLhQM2bMCGFk4Y3z0Zmzzz7b9/8DBw7U0KFDlZ6ertdee02nnXZaCCNznylTpujTTz/V+++/3+g5zkMEU7upSPTo0UORkZGNsuxt27Y1ysZxYLGxsRo0aJBKSkpCHUpYarjihfPRrpSUFKWnp3Ne7mPq1Kl65ZVX9O677+rQQw/1rec8RFtoN4lEdHS0hgwZ0ugnu5cuXaoRI0aEKKrwVVVVpXXr1iklJSXUoYSljIwMJScn+52P1dXVWr58OeejA9u3b1dpaSnn5f8xxmjKlClavHix3nnnHWVkZPg9z3mIttCuhjZmzJihiRMnaujQoRo+fLjmzZunzZs3a/LkyaEOzfWuuuoqjR8/Xr1799a2bdt02223qaKiQrm5uaEOzbV27dqljRs3+h5v2rRJRUVFSkxMVO/evTV9+nTNnj1bffv2Vd++fTV79mx16dJF5513Xgijdpf9HcPExETl5+fr9NNPV0pKir7++mtdf/316tGjh37/+9+HMGr3yMvL09NPP62XX35ZcXFxvspDQkKCOnfuLI/Hw3mI4AvpNSNB4PV6TXp6uomOjjbHHHOM7zIo7N/ZZ59tUlJSTFRUlElNTTWnnXaaWbt2bajDcrV3333XSGq05ObmGmP2XHo3a9Ys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",
      "text/plain": [
       "Graphics object consisting of 1 graphics primitive"
      ]
     },
     "execution_count": 20,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#plot the result of U_complex*U_complex.H to see how far it is from being unitary over the complex numbers\n",
    "plot_arg_complex(gram_rounded, title=f\"arg of gram matrix of complexified uDFT of S_{n} over F_{q**2}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "187d0d5b-3f3c-4d75-866d-36e3a509b54d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[6, 5*z2 + 3, 4*z14^13 + 6*z14^12 + 4*z14^11 + 5*z14^10 + 4*z14^8 + 3*z14^6 + 3*z14^5 + 5*z14^4 + z14^3 + 4*z14^2 + 4*z14 + 1, 2*z14^13 + z14^12 + 2*z14^11 + 5*z14^10 + 6*z14^8 + 2*z14^7 + 6*z14^6 + 5*z14^5 + 4*z14^4 + 3*z14^3 + 3*z14^2 + 4*z14 + 2, z14^13 + 2*z14^12 + 2*z14^11 + 4*z14^10 + 4*z14^9 + 2*z14^6 + z14^3 + 2*z14 + 3, z14^13 + 4*z14^12 + 3*z14^11 + 5*z14^10 + z14^9 + 5*z14^6 + 4*z14^5 + 5*z14^4 + 5*z14^3 + z14^2 + 5*z14 + 3, 3*z14^13 + 6*z14^12 + z14^11 + 3*z14^10 + 4*z14^9 + 6*z14^7 + 5*z14^6 + 3*z14^5 + 6*z14^4 + 2*z14^3 + 4*z14 + 5, 2*z14^13 + 2*z14^12 + z14^11 + 3*z14^9 + 3*z14^8 + 3*z14^7 + 5*z14^5 + 4*z14^4 + 2*z14^3 + 4*z14^2 + 4*z14 + 6, 6*z14^13 + 3*z14^12 + 6*z14^11 + 2*z14^9 + 4*z14^8 + 5*z14^6 + 4*z14^5 + z14^4 + 3*z14^3 + 3*z14^2, 6*z30^29 + 3*z30^28 + 6*z30^27 + 2*z30^25 + 3*z30^24 + 3*z30^23 + 5*z30^22 + 3*z30^21 + 4*z30^20 + 4*z30^19 + 2*z30^17 + 6*z30^16 + 4*z30^15 + 6*z30^14 + 4*z30^12 + 3*z30^11 + 4*z30^10 + 6*z30^9 + 6*z30^8 + 5*z30^7 + 5*z30^6 + 4*z30^5 + 5*z30^4 + 2*z30^3 + 3*z30^2 + 1, z30^29 + 6*z30^28 + 4*z30^26 + 2*z30^25 + 6*z30^24 + 6*z30^23 + 6*z30^22 + 3*z30^21 + z30^20 + 4*z30^19 + 6*z30^18 + 2*z30^17 + z30^16 + 3*z30^15 + 5*z30^13 + 2*z30^12 + 5*z30^11 + 6*z30^10 + 6*z30^9 + 2*z30^8 + 2*z30^6 + 5*z30^5 + 2*z30^4 + 6*z30^3 + 3*z30^2 + z30 + 2, 2*z30^29 + 4*z30^28 + 5*z30^27 + z30^26 + 3*z30^25 + 5*z30^24 + 5*z30^23 + 6*z30^22 + 2*z30^21 + z30^20 + 4*z30^19 + 3*z30^17 + 3*z30^16 + 4*z30^14 + 2*z30^13 + z30^12 + 3*z30^10 + 5*z30^9 + 2*z30^8 + 2*z30^7 + 6*z30^6 + 5*z30^5 + 3*z30^4 + 4*z30^3 + 4*z30^2 + 3*z30 + 2, 4*z30^29 + 5*z30^28 + 6*z30^27 + 2*z30^26 + 2*z30^25 + 5*z30^24 + z30^23 + 6*z30^22 + 6*z30^20 + 3*z30^19 + z30^18 + 4*z30^17 + 4*z30^16 + z30^15 + z30^14 + 6*z30^13 + 3*z30^12 + 3*z30^11 + 2*z30^10 + z30^9 + 2*z30^8 + 5*z30^7 + 2*z30^6 + 6*z30^5 + 4*z30^3 + z30^2 + 4*z30 + 2, 6*z30^29 + 6*z30^28 + 6*z30^27 + z30^26 + 5*z30^25 + 4*z30^24 + 2*z30^23 + 5*z30^22 + 3*z30^21 + 6*z30^20 + z30^18 + 2*z30^17 + 3*z30^16 + 4*z30^15 + 6*z30^14 + 3*z30^13 + 4*z30^12 + 2*z30^11 + 5*z30^10 + z30^8 + 2*z30^7 + 5*z30^6 + 6*z30^4 + 6*z30^3 + 4*z30 + 2, z30^29 + 6*z30^28 + 2*z30^27 + 5*z30^26 + 6*z30^25 + 2*z30^24 + 3*z30^23 + 5*z30^21 + 4*z30^19 + 6*z30^18 + z30^17 + 4*z30^16 + 2*z30^15 + 4*z30^14 + z30^13 + 5*z30^12 + 4*z30^11 + 4*z30^10 + 6*z30^9 + 6*z30^8 + 6*z30^7 + 2*z30^6 + 3*z30^4 + 5*z30^3 + z30^2 + z30 + 3, 4*z30^29 + 6*z30^27 + 6*z30^26 + 4*z30^25 + z30^24 + z30^22 + 6*z30^21 + 6*z30^20 + 5*z30^19 + 5*z30^17 + 4*z30^16 + 3*z30^15 + 3*z30^14 + 6*z30^13 + 6*z30^12 + z30^11 + 2*z30^10 + 4*z30^9 + 5*z30^7 + 5*z30^6 + 3*z30^4 + 5*z30^3 + 2*z30 + 3, 4*z30^29 + z30^28 + 2*z30^27 + 6*z30^24 + 2*z30^22 + 4*z30^21 + 5*z30^20 + 5*z30^19 + 2*z30^17 + z30^16 + z30^13 + 4*z30^12 + 6*z30^11 + 2*z30^10 + 5*z30^9 + 2*z30^8 + 2*z30^7 + z30^6 + z30^5 + 6*z30^4 + z30^3 + 2*z30^2 + 4, 6*z30^29 + 6*z30^28 + 4*z30^27 + 6*z30^26 + z30^25 + z30^24 + 5*z30^23 + 5*z30^22 + 4*z30^21 + 3*z30^20 + 6*z30^19 + 5*z30^18 + 6*z30^17 + z30^16 + z30^15 + 5*z30^14 + 5*z30^13 + 4*z30^12 + 2*z30^11 + 3*z30^10 + 2*z30^8 + 5*z30^7 + z30^6 + 3*z30^4 + 5*z30^3 + 4, 4*z30^29 + 4*z30^28 + 6*z30^27 + 3*z30^26 + 4*z30^25 + z30^24 + 6*z30^22 + 3*z30^21 + 6*z30^20 + 4*z30^19 + z30^18 + 5*z30^17 + 4*z30^16 + 3*z30^15 + 4*z30^13 + z30^12 + 3*z30^11 + 3*z30^10 + 3*z30^9 + z30^8 + 5*z30^7 + z30^5 + 6*z30^4 + 5*z30^3 + z30^2 + 3*z30 + 5, z30^29 + 3*z30^27 + 2*z30^26 + 2*z30^24 + 2*z30^23 + 5*z30^22 + 5*z30^21 + 5*z30^20 + 3*z30^19 + 5*z30^18 + 3*z30^16 + 2*z30^14 + 5*z30^13 + 3*z30^12 + 5*z30^11 + 4*z30^10 + 4*z30^9 + 6*z30^8 + 6*z30^7 + 5*z30^6 + 6*z30^5 + z30^3 + 5*z30^2 + 4*z30 + 5, z30^29 + 4*z30^27 + 2*z30^26 + 4*z30^25 + z30^24 + 5*z30^23 + 6*z30^21 + 4*z30^20 + 4*z30^19 + 3*z30^18 + 4*z30^17 + z30^16 + z30^14 + 3*z30^13 + 4*z30^12 + z30^11 + 2*z30^9 + z30^8 + 4*z30^7 + z30^6 + z30^3 + 5*z30 + 6, 6*z30^28 + 6*z30^27 + 5*z30^26 + 2*z30^25 + 3*z30^23 + 4*z30^21 + z30^20 + 3*z30^19 + 2*z30^18 + 6*z30^17 + 3*z30^16 + 2*z30^15 + 5*z30^14 + 3*z30^13 + 4*z30^9 + 6*z30^8 + 2*z30^7 + 4*z30^6 + 3*z30^5 + z30^4 + z30^2 + 6*z30 + 1, 5*z30^28 + 2*z30^24 + z30^21 + 2*z30^20 + 4*z30^19 + 3*z30^18 + 4*z30^17 + z30^16 + z30^15 + 5*z30^14 + 4*z30^12 + 5*z30^11 + 4*z30^9 + 6*z30^8 + 3*z30^7 + 4*z30^5 + 4*z30^4 + 6*z30^3 + 2*z30 + 2, 4*z30^28 + z30^27 + 2*z30^26 + 5*z30^25 + z30^24 + 4*z30^23 + 6*z30^22 + 3*z30^21 + 3*z30^20 + 5*z30^19 + z30^18 + 6*z30^15 + z30^14 + z30^12 + 6*z30^11 + z30^10 + 4*z30^9 + 5*z30^8 + z30^7 + 2*z30^6 + 3*z30^5 + 2*z30^4 + 3*z30^3 + 5*z30^2 + 4*z30 + 4]"
      ]
     },
     "execution_count": 21,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#compute the eigenvalues of the uDFT matrix\n",
    "eigenvalues = U.eigenvalues(); eigenvalues"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "c61d4b5b-06fb-459a-af3f-a2772a72ad59",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "30"
      ]
     },
     "execution_count": 22,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#compute the maximum degree required for the eigenvalues\n",
    "#this should agree with the splitting field degree, but sometimes it doesn't\n",
    "max_deg_eigs = max([eig.minpoly().degree() for eig in eigenvalues]); max_deg_eigs"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "id": "787c4350-2f3a-4e24-8e48-b664ea502f44",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "True"
      ]
     },
     "execution_count": 23,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#compute the multiplicity of each eigenvalue. check if all eigenvalues are distinct\n",
    "from collections import Counter\n",
    "multiplicities = Counter(eigenvalues)\n",
    "all_unique = all(count == 1 for count in multiplicities.values()); all_unique"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 24,
   "id": "577cf124-887a-452b-ad8e-d218ab5a8cbb",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "x^24 + z2*x^23 + (2*z2 + 2)*x^22 + (4*z2 + 4)*x^21 + (2*z2 + 6)*x^20 + 4*z2*x^19 + 2*z2*x^18 + (4*z2 + 3)*x^17 + (4*z2 + 6)*x^16 + (5*z2 + 2)*x^15 + 2*z2*x^14 + (6*z2 + 2)*x^13 + 5*z2*x^11 + (z2 + 5)*x^10 + (z2 + 1)*x^9 + z2*x^8 + (5*z2 + 5)*x^7 + (z2 + 5)*x^6 + (2*z2 + 3)*x^5 + 2*x^4 + (6*z2 + 1)*x^3 + (3*z2 + 4)*x^2 + (4*z2 + 6)*x + z2 + 3"
      ]
     },
     "execution_count": 24,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#compute the charpoly of the uDFT matrix \n",
    "charpoly = U.minimal_polynomial(); charpoly"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "id": "9ec1d71f-fd57-4261-9129-a634f0184e24",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "field containing all eigenvalues: K = Finite Field in z30 of size 7^30\n",
      "splitting field: L = Finite Field in a of size 7^210\n",
      "minpoly of K:  x^30 + x^22 + 4*x^21 + 4*x^20 + x^19 + 4*x^18 + x^16 + 2*x^15 + 3*x^14 + 6*x^13 + 5*x^12 + 2*x^11 + 3*x^10 + 3*x^9 + 2*x^8 + 4*x^7 + 2*x^6 + 3*x^5 + x^3 + 5*x^2 + 2*x + 3\n"
     ]
    }
   ],
   "source": [
    "#compute a splitting field of the characteristic polynomial\n",
    "K = GF(q**max_deg_eigs); print(f\"field containing all eigenvalues: K = {K}\")\n",
    "L = charpoly.splitting_field('a'); print(f\"splitting field: L = {L}\")\n",
    "print(\"minpoly of K: \", K.multiplicative_generator().minimal_polynomial())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 26,
   "id": "6ad9735a-b9fb-46a5-822e-088f724c87a5",
   "metadata": {},
   "outputs": [],
   "source": [
    "def save_array(array, path, filename):\n",
    "    \"\"\"\n",
    "    save an array as a comma separated value file by converting the elements to strings\n",
    "    \"\"\"\n",
    "    import os\n",
    "    if not os.path.exists(path):\n",
    "        os.makedirs(path)  # Create the directory if it doesn't exist\n",
    "    full_path = os.path.join(path, filename)\n",
    "    with open(full_path, 'w') as f:\n",
    "        for i, element in enumerate(array):\n",
    "            f.write(str(element))\n",
    "            if i < len(array) - 1:\n",
    "                f.write(\",\")\n",
    "        f.write(\"\\n\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "id": "a8d24cd3-3d8e-4e60-ab01-63a3fac92936",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute and save the eigenvalues over K\n",
    "eigenvalues_K = matrix(K,U).eigenvalues(extend=False)\n",
    "save_array(eigenvalues_K, path='data/eigenvalues', filename=f\"eigenvalues_uDFT_n={n}_q={q}_deg={K.degree()}\" + '.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 37,
   "id": "2ce01f0b-1c29-4761-a1a5-8ec4b30aeb90",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute and save the eigenvalues over L\n",
    "eigenvalues_L = matrix(L,U).eigenvalues(extend=False)\n",
    "save_array(eigenvalues_L, path='data/eigenvalues', filename=f\"eigenvalues_uDFT_n={n}_q={q}_deg={L.degree()}\" + '.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 38,
   "id": "d1b9c85e-3854-4210-90d4-6dd27c07e836",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute and save the discrete log of the eigenvalues over K\n",
    "log_eigenvalues_K = list(map(lambda x: discrete_log(K,x), eigenvalues_K))\n",
    "save_array(log_eigenvalues_K, path='data/discrete_log_eigenvalues', filename=f\"discrete_log_eigenvalues_uDFT_n={n}_q={q}_deg={K.degree()}\" + '.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 127,
   "id": "a43cd6e7-9147-4709-bd33-213d71716b28",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute and save the discrete log of the eigenvalues over L\n",
    "log_eigenvalues_L = list(map(lambda x: discrete_log(L,x), eigenvalues_L))\n",
    "save_array(log_eigenvalues_L, path='data/discrete_log_eigenvalues', filename=f\"discrete_log_eigenvalues_uDFT_n={n}_q={q}_deg={L.degree()}\" + '.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 29,
   "id": "c41fcb01-39a8-4e24-a828-45a0f19bf842",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute and save the complexified eigenvalue over K\n",
    "complexified_eigenvalues_K = [brauer_map(K,log_a=log_eig) for log_eig in log_eigenvalues_K]\n",
    "save_array(complexified_eigenvalues_K, path='data/complexified_eigenvalues', filename=f\"complexified_eigenvalues_uDFT_n={n}_q={q}_deg={K.degree()}\" + '.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 32,
   "id": "271cd002-0f0e-4b8f-beb4-8e3a02d0c959",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute and save the complexified eigenvalues over L\n",
    "complexified_eigenvalues_L = [brauer_map(L,log_a=log_eig) for log_eig in log_eigenvalues_L]\n",
    "save_array(complexified_eigenvalues_L, path='data/complexified_eigenvalues', filename=f\"complexified_eigenvalues_uDFT_n={n}_q={q}_deg={L.degree()}\" + '.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 30,
   "id": "997124ab-f397-42cf-a36d-fc1670e27bf0",
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "import numpy as np\n",
    "\n",
    "# Extract real and imaginary parts for plotting\n",
    "real_parts = [eig.real() for eig in complexified_eigenvalues_K]\n",
    "imaginary_parts = [eig.imag() for eig in complexified_eigenvalues_K]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 36,
   "id": "92b0d6ae-51a8-43d2-9ead-f4b40d89b9ac",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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",
      "text/plain": [
       "<Figure size 600x600 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# Plotting the eigenvalues on the complex plane\n",
    "title = f\"complexified eigenvalues of uDFT for n={n} q={q} deg={K.degree()}\"\n",
    "plt.figure(figsize=(6,6))\n",
    "plt.scatter(real_parts, imaginary_parts, color='blue', label=\"Eigenvalues\")\n",
    "plt.axhline(0, color='black',linewidth=0.5)\n",
    "plt.axvline(0, color='black',linewidth=0.5)\n",
    "plt.gca().set_aspect('equal', adjustable='box')\n",
    "plt.title(title)\n",
    "plt.xlabel('Real Part')\n",
    "plt.ylabel('Imaginary Part')\n",
    "plt.grid(True)\n",
    "plt.savefig('plots/eigenvalues/finite_fields/' + title.replace(' ','_') + '.png', dpi=300, bbox_inches=\"tight\")\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 32,
   "id": "f4558cc4-ca09-4785-affd-9da688591bb2",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute the eigenvectors over a splitting field\n",
    "D, P = matrix(L,U).eigenmatrix_right()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 33,
   "id": "25a882fa-3dd4-4997-857d-bc6b553fc66a",
   "metadata": {},
   "outputs": [],
   "source": [
    "def write_matrix_to_csv(P,filename):\n",
    "    \"\"\"\n",
    "    dump the eigenvector matrix into a .csv file\n",
    "    \"\"\"\n",
    "    with open(filename, \"w\") as f:\n",
    "        for row in P:\n",
    "            f.write(\",\".join(map(str, row)) + \"\\n\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 34,
   "id": "852531ea-24de-4d96-8bd6-192b399c04b8",
   "metadata": {},
   "outputs": [],
   "source": [
    "write_matrix_to_csv(P,filename=f\"data/eigenvectors/eigenvectors_unitary_dft_symmetric_group_finite_field_n={n}_q={q}.csv\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 323,
   "id": "c302f2a7-4bda-4203-8ab5-dca8f1f173d4",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[                0                 0                 0                 0                 0                 0]\n",
       "[22689761797238269 20702645593902589 12772541917661592 13523245755449361  3495018369229650 11976932273987601]\n",
       "[11727754551588639  6978876427228959 21002152480419206  1362287283069904 15409676675944734 20973400951923664]\n",
       "[14579426119294872  7816719458901912  1589044626709508 20061842302069898  3404489112688662  3103075099765898]\n",
       "[ 5077809292629124  7357597852275844 17981569188489612 11347564654638395  2911937431490442 19959385578951995]\n",
       "[16327908218895593  5293676504325353 12706162947981234 10976243789732849  6310175695144002  3732098528082929]"
      ]
     },
     "execution_count": 323,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#take the discrete log of the eigenvector matrix\n",
    "log_P = P.apply_map(lambda x: discrete_log(L,x)); log_P"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 46,
   "id": "0b3f2964-f2a2-44bb-ab85-70a176f10929",
   "metadata": {},
   "outputs": [
    {
     "ename": "TypeError",
     "evalue": "Image data of dtype object cannot be converted to float",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mTypeError\u001b[0m                                 Traceback (most recent call last)",
      "Cell \u001b[0;32mIn[46], line 2\u001b[0m\n\u001b[1;32m      1\u001b[0m \u001b[38;5;66;03m#plot the discrete log of the eigenvector matrix P\u001b[39;00m\n\u001b[0;32m----> 2\u001b[0m \u001b[43mplot_discrete_log\u001b[49m\u001b[43m(\u001b[49m\u001b[43mL\u001b[49m\u001b[43m,\u001b[49m\u001b[43mlog_P\u001b[49m\u001b[43m,\u001b[49m\u001b[43mpath\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[38;5;124;43m'\u001b[39;49m\u001b[38;5;124;43mplots/eigenvectors/\u001b[39;49m\u001b[38;5;124;43m'\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43mtitle\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[38;5;124;43mf\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[38;5;124;43mdiscrete log of eigenvectors of uDFT for n=\u001b[39;49m\u001b[38;5;132;43;01m{\u001b[39;49;00m\u001b[43mn\u001b[49m\u001b[38;5;132;43;01m}\u001b[39;49;00m\u001b[38;5;124;43m q=\u001b[39;49m\u001b[38;5;132;43;01m{\u001b[39;49;00m\u001b[43mq\u001b[49m\u001b[38;5;132;43;01m}\u001b[39;49;00m\u001b[38;5;124;43m deg=\u001b[39;49m\u001b[38;5;132;43;01m{\u001b[39;49;00m\u001b[43mL\u001b[49m\u001b[38;5;241;43m.\u001b[39;49m\u001b[43mdegree\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\u001b[38;5;132;43;01m}\u001b[39;49;00m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\n",
      "Cell \u001b[0;32mIn[10], line 22\u001b[0m, in \u001b[0;36mplot_discrete_log\u001b[0;34m(F, M, path, title)\u001b[0m\n\u001b[1;32m     19\u001b[0m norm \u001b[38;5;241m=\u001b[39m BoundaryNorm([\u001b[38;5;241m-\u001b[39mInteger(\u001b[38;5;241m1\u001b[39m)] \u001b[38;5;241m+\u001b[39m \u001b[38;5;28mlist\u001b[39m(np\u001b[38;5;241m.\u001b[39mlinspace(Integer(\u001b[38;5;241m0\u001b[39m), F\u001b[38;5;241m.\u001b[39morder()\u001b[38;5;241m-\u001b[39mInteger(\u001b[38;5;241m1\u001b[39m), num_colors)), custom_cmap\u001b[38;5;241m.\u001b[39mN) \u001b[38;5;66;03m#map -1 to black\u001b[39;00m\n\u001b[1;32m     21\u001b[0m \u001b[38;5;66;03m# Plotting the data\u001b[39;00m\n\u001b[0;32m---> 22\u001b[0m \u001b[43mplt\u001b[49m\u001b[38;5;241;43m.\u001b[39;49m\u001b[43mimshow\u001b[49m\u001b[43m(\u001b[49m\u001b[43mM\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mcmap\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mcustom_cmap\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mnorm\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mnorm\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43minterpolation\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[38;5;124;43mnearest\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[43m)\u001b[49m\n\u001b[1;32m     23\u001b[0m plt\u001b[38;5;241m.\u001b[39mtitle(title, fontsize\u001b[38;5;241m=\u001b[39mInteger(\u001b[38;5;241m16\u001b[39m))\n\u001b[1;32m     24\u001b[0m plt\u001b[38;5;241m.\u001b[39mcolorbar()\n",
      "File \u001b[0;32m/private/var/tmp/sage-10.5-current/local/var/lib/sage/venv-python3.12.5/lib/python3.12/site-packages/matplotlib/pyplot.py:3346\u001b[0m, in \u001b[0;36mimshow\u001b[0;34m(X, cmap, norm, aspect, interpolation, alpha, vmin, vmax, origin, extent, interpolation_stage, filternorm, filterrad, resample, url, data, **kwargs)\u001b[0m\n\u001b[1;32m   3325\u001b[0m \u001b[38;5;129m@_copy_docstring_and_deprecators\u001b[39m(Axes\u001b[38;5;241m.\u001b[39mimshow)\n\u001b[1;32m   3326\u001b[0m \u001b[38;5;28;01mdef\u001b[39;00m \u001b[38;5;21mimshow\u001b[39m(\n\u001b[1;32m   3327\u001b[0m     X: ArrayLike \u001b[38;5;241m|\u001b[39m PIL\u001b[38;5;241m.\u001b[39mImage\u001b[38;5;241m.\u001b[39mImage,\n\u001b[0;32m   (...)\u001b[0m\n\u001b[1;32m   3344\u001b[0m     \u001b[38;5;241m*\u001b[39m\u001b[38;5;241m*\u001b[39mkwargs,\n\u001b[1;32m   3345\u001b[0m ) \u001b[38;5;241m-\u001b[39m\u001b[38;5;241m>\u001b[39m AxesImage:\n\u001b[0;32m-> 3346\u001b[0m     __ret \u001b[38;5;241m=\u001b[39m \u001b[43mgca\u001b[49m\u001b[43m(\u001b[49m\u001b[43m)\u001b[49m\u001b[38;5;241;43m.\u001b[39;49m\u001b[43mimshow\u001b[49m\u001b[43m(\u001b[49m\n\u001b[1;32m   3347\u001b[0m \u001b[43m        \u001b[49m\u001b[43mX\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3348\u001b[0m \u001b[43m        \u001b[49m\u001b[43mcmap\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mcmap\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3349\u001b[0m \u001b[43m        \u001b[49m\u001b[43mnorm\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mnorm\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3350\u001b[0m \u001b[43m        \u001b[49m\u001b[43maspect\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43maspect\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3351\u001b[0m \u001b[43m        \u001b[49m\u001b[43minterpolation\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43minterpolation\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3352\u001b[0m \u001b[43m        \u001b[49m\u001b[43malpha\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43malpha\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3353\u001b[0m \u001b[43m        \u001b[49m\u001b[43mvmin\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mvmin\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3354\u001b[0m \u001b[43m        \u001b[49m\u001b[43mvmax\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mvmax\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3355\u001b[0m \u001b[43m        \u001b[49m\u001b[43morigin\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43morigin\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3356\u001b[0m \u001b[43m        \u001b[49m\u001b[43mextent\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mextent\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3357\u001b[0m \u001b[43m        \u001b[49m\u001b[43minterpolation_stage\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43minterpolation_stage\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3358\u001b[0m \u001b[43m        \u001b[49m\u001b[43mfilternorm\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mfilternorm\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3359\u001b[0m \u001b[43m        \u001b[49m\u001b[43mfilterrad\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mfilterrad\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3360\u001b[0m \u001b[43m        \u001b[49m\u001b[43mresample\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mresample\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3361\u001b[0m \u001b[43m        \u001b[49m\u001b[43murl\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43murl\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3362\u001b[0m \u001b[43m        \u001b[49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[43m(\u001b[49m\u001b[43m{\u001b[49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[38;5;124;43mdata\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[43m:\u001b[49m\u001b[43m \u001b[49m\u001b[43mdata\u001b[49m\u001b[43m}\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;28;43;01mif\u001b[39;49;00m\u001b[43m \u001b[49m\u001b[43mdata\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;129;43;01mis\u001b[39;49;00m\u001b[43m \u001b[49m\u001b[38;5;129;43;01mnot\u001b[39;49;00m\u001b[43m \u001b[49m\u001b[38;5;28;43;01mNone\u001b[39;49;00m\u001b[43m \u001b[49m\u001b[38;5;28;43;01melse\u001b[39;49;00m\u001b[43m \u001b[49m\u001b[43m{\u001b[49m\u001b[43m}\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3363\u001b[0m \u001b[43m        \u001b[49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[43mkwargs\u001b[49m\u001b[43m,\u001b[49m\n\u001b[1;32m   3364\u001b[0m \u001b[43m    \u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m   3365\u001b[0m     sci(__ret)\n\u001b[1;32m   3366\u001b[0m     \u001b[38;5;28;01mreturn\u001b[39;00m __ret\n",
      "File \u001b[0;32m/private/var/tmp/sage-10.5-current/local/var/lib/sage/venv-python3.12.5/lib/python3.12/site-packages/matplotlib/__init__.py:1465\u001b[0m, in \u001b[0;36m_preprocess_data.<locals>.inner\u001b[0;34m(ax, data, *args, **kwargs)\u001b[0m\n\u001b[1;32m   1462\u001b[0m \u001b[38;5;129m@functools\u001b[39m\u001b[38;5;241m.\u001b[39mwraps(func)\n\u001b[1;32m   1463\u001b[0m \u001b[38;5;28;01mdef\u001b[39;00m \u001b[38;5;21minner\u001b[39m(ax, \u001b[38;5;241m*\u001b[39margs, data\u001b[38;5;241m=\u001b[39m\u001b[38;5;28;01mNone\u001b[39;00m, \u001b[38;5;241m*\u001b[39m\u001b[38;5;241m*\u001b[39mkwargs):\n\u001b[1;32m   1464\u001b[0m     \u001b[38;5;28;01mif\u001b[39;00m data \u001b[38;5;129;01mis\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m:\n\u001b[0;32m-> 1465\u001b[0m         \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mfunc\u001b[49m\u001b[43m(\u001b[49m\u001b[43max\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[38;5;28;43mmap\u001b[39;49m\u001b[43m(\u001b[49m\u001b[43msanitize_sequence\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43margs\u001b[49m\u001b[43m)\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[43mkwargs\u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m   1467\u001b[0m     bound \u001b[38;5;241m=\u001b[39m new_sig\u001b[38;5;241m.\u001b[39mbind(ax, \u001b[38;5;241m*\u001b[39margs, \u001b[38;5;241m*\u001b[39m\u001b[38;5;241m*\u001b[39mkwargs)\n\u001b[1;32m   1468\u001b[0m     auto_label \u001b[38;5;241m=\u001b[39m (bound\u001b[38;5;241m.\u001b[39marguments\u001b[38;5;241m.\u001b[39mget(label_namer)\n\u001b[1;32m   1469\u001b[0m                   \u001b[38;5;129;01mor\u001b[39;00m bound\u001b[38;5;241m.\u001b[39mkwargs\u001b[38;5;241m.\u001b[39mget(label_namer))\n",
      "File \u001b[0;32m/private/var/tmp/sage-10.5-current/local/var/lib/sage/venv-python3.12.5/lib/python3.12/site-packages/matplotlib/axes/_axes.py:5751\u001b[0m, in \u001b[0;36mAxes.imshow\u001b[0;34m(self, X, cmap, norm, aspect, interpolation, alpha, vmin, vmax, origin, extent, interpolation_stage, filternorm, filterrad, resample, url, **kwargs)\u001b[0m\n\u001b[1;32m   5748\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m aspect \u001b[38;5;129;01mis\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m:\n\u001b[1;32m   5749\u001b[0m     \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39mset_aspect(aspect)\n\u001b[0;32m-> 5751\u001b[0m \u001b[43mim\u001b[49m\u001b[38;5;241;43m.\u001b[39;49m\u001b[43mset_data\u001b[49m\u001b[43m(\u001b[49m\u001b[43mX\u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m   5752\u001b[0m im\u001b[38;5;241m.\u001b[39mset_alpha(alpha)\n\u001b[1;32m   5753\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m im\u001b[38;5;241m.\u001b[39mget_clip_path() \u001b[38;5;129;01mis\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m:\n\u001b[1;32m   5754\u001b[0m     \u001b[38;5;66;03m# image does not already have clipping set, clip to axes patch\u001b[39;00m\n",
      "File \u001b[0;32m/private/var/tmp/sage-10.5-current/local/var/lib/sage/venv-python3.12.5/lib/python3.12/site-packages/matplotlib/image.py:723\u001b[0m, in \u001b[0;36m_ImageBase.set_data\u001b[0;34m(self, A)\u001b[0m\n\u001b[1;32m    721\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;28misinstance\u001b[39m(A, PIL\u001b[38;5;241m.\u001b[39mImage\u001b[38;5;241m.\u001b[39mImage):\n\u001b[1;32m    722\u001b[0m     A \u001b[38;5;241m=\u001b[39m pil_to_array(A)  \u001b[38;5;66;03m# Needed e.g. to apply png palette.\u001b[39;00m\n\u001b[0;32m--> 723\u001b[0m \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39m_A \u001b[38;5;241m=\u001b[39m \u001b[38;5;28;43mself\u001b[39;49m\u001b[38;5;241;43m.\u001b[39;49m\u001b[43m_normalize_image_array\u001b[49m\u001b[43m(\u001b[49m\u001b[43mA\u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m    724\u001b[0m \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39m_imcache \u001b[38;5;241m=\u001b[39m \u001b[38;5;28;01mNone\u001b[39;00m\n\u001b[1;32m    725\u001b[0m \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39mstale \u001b[38;5;241m=\u001b[39m \u001b[38;5;28;01mTrue\u001b[39;00m\n",
      "File \u001b[0;32m/private/var/tmp/sage-10.5-current/local/var/lib/sage/venv-python3.12.5/lib/python3.12/site-packages/matplotlib/image.py:688\u001b[0m, in \u001b[0;36m_ImageBase._normalize_image_array\u001b[0;34m(A)\u001b[0m\n\u001b[1;32m    686\u001b[0m A \u001b[38;5;241m=\u001b[39m cbook\u001b[38;5;241m.\u001b[39msafe_masked_invalid(A, copy\u001b[38;5;241m=\u001b[39m\u001b[38;5;28;01mTrue\u001b[39;00m)\n\u001b[1;32m    687\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m A\u001b[38;5;241m.\u001b[39mdtype \u001b[38;5;241m!=\u001b[39m np\u001b[38;5;241m.\u001b[39muint8 \u001b[38;5;129;01mand\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m np\u001b[38;5;241m.\u001b[39mcan_cast(A\u001b[38;5;241m.\u001b[39mdtype, \u001b[38;5;28mfloat\u001b[39m, \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124msame_kind\u001b[39m\u001b[38;5;124m\"\u001b[39m):\n\u001b[0;32m--> 688\u001b[0m     \u001b[38;5;28;01mraise\u001b[39;00m \u001b[38;5;167;01mTypeError\u001b[39;00m(\u001b[38;5;124mf\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mImage data of dtype \u001b[39m\u001b[38;5;132;01m{\u001b[39;00mA\u001b[38;5;241m.\u001b[39mdtype\u001b[38;5;132;01m}\u001b[39;00m\u001b[38;5;124m cannot be \u001b[39m\u001b[38;5;124m\"\u001b[39m\n\u001b[1;32m    689\u001b[0m                     \u001b[38;5;124mf\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mconverted to float\u001b[39m\u001b[38;5;124m\"\u001b[39m)\n\u001b[1;32m    690\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m A\u001b[38;5;241m.\u001b[39mndim \u001b[38;5;241m==\u001b[39m \u001b[38;5;241m3\u001b[39m \u001b[38;5;129;01mand\u001b[39;00m A\u001b[38;5;241m.\u001b[39mshape[\u001b[38;5;241m-\u001b[39m\u001b[38;5;241m1\u001b[39m] \u001b[38;5;241m==\u001b[39m \u001b[38;5;241m1\u001b[39m:\n\u001b[1;32m    691\u001b[0m     A \u001b[38;5;241m=\u001b[39m A\u001b[38;5;241m.\u001b[39msqueeze(\u001b[38;5;241m-\u001b[39m\u001b[38;5;241m1\u001b[39m)  \u001b[38;5;66;03m# If just (M, N, 1), assume scalar and apply colormap.\u001b[39;00m\n",
      "\u001b[0;31mTypeError\u001b[0m: Image data of dtype object cannot be converted to float"
     ]
    },
    {
     "data": {
      "image/png": "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5z//WXfccYckaeXKlfr000/1/vvvR/bft2+f1qxZo48++kiTJ0/WunXrIsEDACAWMQUMAIDRxu9CBACYRMAAACYRMACASQQMAGDSZRMw/kTL4HhZp127dmn+/Pm6+uqrlZSUpLlz5+q9994bwdmODq8/S30OHjwov9+v2267bXgneJnwuk7hcFilpaVKS0tTIBDQDTfcoMrKyhGa7ejxuk5VVVWaNWuWrrjiCqWkpOiBBx5QZ2fnCM12dOzfv1+LFy/W5MmT5fP59M477/zgMUPyGj6kv1kxRn/961/duHHj3Ouvv+4aGxvdY4895iZMmOA+++yzAfc/duyYu+KKK9xjjz3mGhsb3euvv+7GjRvn3nrrrRGe+cjyuk6PPfaYe+aZZ9y///1v9/HHH7uSkhI3btw499///neEZz5yvK5Rn1OnTrlp06a5vLw8N2vWrJGZ7CiKZZ3uu+8+l52d7UKhkDt+/Lj717/+5Q4ePDiCsx55XteptrbWxcXFueeff94dO3bM1dbWultuucUtWbJkhGc+smpqalxpaal7++23nSS3e/fu791/qF7DL4uAzZkzxxUWFkaN3XTTTW79+vUD7v+HP/zB3XTTTVFjDz30kLv99tuHbY6XA6/rNJCbb77Zbdy4caindtmIdY3y8/PdH//4R7dhw4YxETCv6/T3v//dBYNB19nZORLTu2x4Xac//elPbtq0aVFjL7zwgps6deqwzfFyM5iADdVr+KhfQjx79qzq6uqUl5cXNZ6Xl6dDhw4NeMzhw4f77X/vvffqyJEj+uabb4ZtrqMplnX6rvPnz6unp2fIfyP05SLWNdq2bZuOHj2qDRs2DPcULwuxrNOePXuUlZWlZ599VlOmTNGMGTO0du1aff311yMx5VERyzrl5OToxIkTqqmpkXNOX375pd566y0tWrRoJKZsxlC9hsf02+iH0kj9iRbrYlmn73ruuef01VdfaenSpcMxxVEXyxp98sknWr9+vWpra+X3j/r/DiMilnU6duyYDhw4oMTERO3evVsdHR16+OGHdfLkyR/t+2CxrFNOTo6qqqqUn5+v//3vfzp37pzuu+8+vfjiiyMxZTOG6jV81M/A+gz3n2j5sfC6Tn127typp556StXV1brmmmuGa3qXhcGuUW9vr5YtW6aNGzdqxowZIzW9y4aXn6Xz58/L5/OpqqpKc+bM0cKFC7V582Zt3779R30WJnlbp8bGRhUVFenJJ59UXV2d3n33XR0/fpxfnTeAoXgNH/V/cvInWgYnlnXqU11drVWrVunNN9/UPffcM5zTHFVe16inp0dHjhxRfX29Hn30UUkXXqidc/L7/dq7d6/uvvvuEZn7SIrlZyklJUVTpkyJ+oOyGRkZcs7pxIkTmj59+rDOeTTEsk5lZWWaN2+eHn/8cUnSrbfeqgkTJig3N1dPP/30j/LqUCyG6jV81M/A+BMtgxPLOkkXzrxWrlypHTt2/Oivw3tdo6SkJH344YdqaGiIbIWFhbrxxhvV0NCg7OzskZr6iIrlZ2nevHn64osvdPr06cjYxx9/rLi4OE2dOnVY5ztaYlmnM2fO9PujjfHx8ZL+/wwDQ/ga7ukjH8Ok76OqW7dudY2NjW716tVuwoQJ7tNPP3XOObd+/Xq3fPnyyP59H8Fcs2aNa2xsdFu3bh1TH6Mf7Drt2LHD+f1+9/LLL7vW1tbIdurUqdF6CsPO6xp911j5FKLXderp6XFTp051v/71r91HH33k9u3b56ZPn+4efPDB0XoKI8LrOm3bts35/X5XXl7ujh496g4cOOCysrLcnDlzRuspjIienh5XX1/v6uvrnSS3efNmV19fH/m6wXC9hl8WAXPOuZdfftmlpaW5hIQEN3v2bLdv377If1uxYoW78847o/Z///333c9//nOXkJDgrr/+eldRUTHCMx4dXtbpzjvvdJL6bStWrBj5iY8grz9L3zZWAuac93Vqampy99xzjxs/frybOnWqKy4udmfOnBnhWY88r+v0wgsvuJtvvtmNHz/epaSkuN/85jfuxIkTIzzrkfWPf/zje19rhus1nD+nAgAwadTfAwMAIBYEDABgEgEDAJhEwAAAJhEwAIBJBAwAYBIBAwCYRMAAACYRMACASQQMAGASAQMAmPR/vVBObw9VdzEAAAAASUVORK5CYII=",
      "text/plain": [
       "<Figure size 640x480 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "#plot the discrete log of the eigenvector matrix P\n",
    "plot_discrete_log(L,log_P,path='plots/eigenvectors/',title=f\"discrete log of eigenvectors of uDFT for n={n} q={q} deg={L.degree()}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 35,
   "id": "67a4f9ac-a80e-4535-bab2-f2c041e16291",
   "metadata": {},
   "outputs": [],
   "source": [
    "def dump_discrete_log(L, P):\n",
    "    \"\"\"\n",
    "    - compute the discrete log of each entry of P, the eigenvector matrix of the uDFT of S_n over F_{q^2}\n",
    "    - write the resulting matrix to a comma separated value file\n",
    "    - include a progress bar since computing discrete logs for large fields takes a long time\n",
    "    \"\"\"\n",
    "    total = P.nrows() * P.ncols()\n",
    "    count = 0\n",
    "    filename = f\"data/eigenvectors/discrete_log_dft_symmetric_group_finite_field_eigenvector_matrix_n={n}_q={q}.csv\"\n",
    "    with open(filename, \"w\") as f:\n",
    "        for i in range(P.nrows()):\n",
    "            for j in range(P.ncols()):\n",
    "                log_value = discrete_log(L, P[i, j])\n",
    "                f.write(str(log_value))\n",
    "                if j != P.ncols()-1:\n",
    "                    f.write(\",\")\n",
    "                count += 1\n",
    "                # Print progress as a percentage\n",
    "                progress = float(count) / total * 100\n",
    "                print(f\"Progress: {progress:.2f}%\\n\", end=\"\\r\", flush=True)\n",
    "            f.write(\"\\n\")\n",
    "    return \"done\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 40,
   "id": "c62f8f4f-2da2-48ee-a4a6-2b93f0d631a9",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Progress: 0.17%\n",
      "Progress: 0.35%\n",
      "Progress: 0.52%\n",
      "Progress: 0.69%\n",
      "Progress: 0.87%\n",
      "Progress: 1.04%\n",
      "Progress: 1.22%\n",
      "Progress: 1.39%\n",
      "Progress: 1.56%\n",
      "Progress: 1.74%\n",
      "Progress: 1.91%\n",
      "Progress: 2.08%\n",
      "Progress: 2.26%\n",
      "Progress: 2.43%\n",
      "Progress: 2.60%\n",
      "Progress: 2.78%\n",
      "Progress: 2.95%\n",
      "Progress: 3.12%\n",
      "Progress: 3.30%\n",
      "Progress: 3.47%\n",
      "Progress: 3.65%\n",
      "Progress: 3.82%\n",
      "Progress: 3.99%\n",
      "Progress: 4.17%\n",
      "Progress: 4.34%\n",
      "Progress: 4.51%\n",
      "Progress: 4.69%\n",
      "Progress: 4.86%\n",
      "Progress: 5.03%\n",
      "Progress: 5.21%\n",
      "Progress: 5.38%\n",
      "Progress: 5.56%\n",
      "Progress: 5.73%\n",
      "Progress: 5.90%\n",
      "Progress: 6.08%\n",
      "Progress: 6.25%\n",
      "Progress: 6.42%\n",
      "Progress: 6.60%\n",
      "Progress: 6.77%\n",
      "Progress: 6.94%\n",
      "Progress: 7.12%\n",
      "Progress: 7.29%\n",
      "Progress: 7.47%\n",
      "Progress: 7.64%\n",
      "Progress: 7.81%\n",
      "Progress: 7.99%\n",
      "Progress: 8.16%\n",
      "Progress: 8.33%\n",
      "Progress: 8.51%\n",
      "Progress: 8.68%\n",
      "Progress: 8.85%\n",
      "Progress: 9.03%\n",
      "Progress: 9.20%\n",
      "Progress: 9.38%\n",
      "Progress: 9.55%\n",
      "Progress: 9.72%\n",
      "Progress: 9.90%\n",
      "Progress: 10.07%\n",
      "Progress: 10.24%\n",
      "Progress: 10.42%\n",
      "Progress: 10.59%\n",
      "Progress: 10.76%\n",
      "Progress: 10.94%\n",
      "Progress: 11.11%\n",
      "Progress: 11.28%\n",
      "Progress: 11.46%\n",
      "Progress: 11.63%\n",
      "Progress: 11.81%\n",
      "Progress: 11.98%\n",
      "Progress: 12.15%\n",
      "Progress: 12.33%\n",
      "Progress: 12.50%\n",
      "Progress: 12.67%\n",
      "Progress: 12.85%\n",
      "Progress: 13.02%\n",
      "Progress: 13.19%\n",
      "Progress: 13.37%\n",
      "Progress: 13.54%\n",
      "Progress: 13.72%\n",
      "Progress: 13.89%\n",
      "Progress: 14.06%\n",
      "Progress: 14.24%\n",
      "Progress: 14.41%\n",
      "Progress: 14.58%\n",
      "Progress: 14.76%\n",
      "Progress: 14.93%\n",
      "Progress: 15.10%\n",
      "Progress: 15.28%\n",
      "Progress: 15.45%\n",
      "Progress: 15.62%\n",
      "Progress: 15.80%\n",
      "Progress: 15.97%\n",
      "Progress: 16.15%\n",
      "Progress: 16.32%\n",
      "Progress: 16.49%\n",
      "Progress: 16.67%\n",
      "Progress: 16.84%\n",
      "Progress: 17.01%\n",
      "Progress: 17.19%\n",
      "Progress: 17.36%\n",
      "Progress: 17.53%\n",
      "Progress: 17.71%\n",
      "Progress: 17.88%\n",
      "Progress: 18.06%\n",
      "Progress: 18.23%\n",
      "Progress: 18.40%\n",
      "Progress: 18.58%\n",
      "Progress: 18.75%\n",
      "Progress: 18.92%\n",
      "Progress: 19.10%\n",
      "Progress: 19.27%\n",
      "Progress: 19.44%\n",
      "Progress: 19.62%\n",
      "Progress: 19.79%\n",
      "Progress: 19.97%\n",
      "Progress: 20.14%\n",
      "Progress: 20.31%\n",
      "Progress: 20.49%\n",
      "Progress: 20.66%\n",
      "Progress: 20.83%\n",
      "Progress: 21.01%\n",
      "Progress: 21.18%\n",
      "Progress: 21.35%\n",
      "Progress: 21.53%\n",
      "Progress: 21.70%\n",
      "Progress: 21.88%\n",
      "Progress: 22.05%\n",
      "Progress: 22.22%\n",
      "Progress: 22.40%\n",
      "Progress: 22.57%\n",
      "Progress: 22.74%\n",
      "Progress: 22.92%\n",
      "Progress: 23.09%\n",
      "Progress: 23.26%\n",
      "Progress: 23.44%\n",
      "Progress: 23.61%\n",
      "Progress: 23.78%\n",
      "Progress: 23.96%\n",
      "Progress: 24.13%\n",
      "Progress: 24.31%\n",
      "Progress: 24.48%\n",
      "Progress: 24.65%\n",
      "Progress: 24.83%\n",
      "Progress: 25.00%\n",
      "Progress: 25.17%\n",
      "Progress: 25.35%\n",
      "Progress: 25.52%\n",
      "Progress: 25.69%\n",
      "Progress: 25.87%\n",
      "Progress: 26.04%\n",
      "Progress: 26.22%\n",
      "Progress: 26.39%\n",
      "Progress: 26.56%\n",
      "Progress: 26.74%\n",
      "Progress: 26.91%\n",
      "Progress: 27.08%\n",
      "Progress: 27.26%\n",
      "Progress: 27.43%\n",
      "Progress: 27.60%\n",
      "Progress: 27.78%\n",
      "Progress: 27.95%\n",
      "Progress: 28.12%\n",
      "Progress: 28.30%\n",
      "Progress: 28.47%\n",
      "Progress: 28.65%\n",
      "Progress: 28.82%\n",
      "Progress: 28.99%\n",
      "Progress: 29.17%\n",
      "Progress: 29.34%\n",
      "Progress: 29.51%\n",
      "Progress: 29.69%\n",
      "Progress: 29.86%\n",
      "Progress: 30.03%\n",
      "Progress: 30.21%\n",
      "Progress: 30.38%\n",
      "Progress: 30.56%\n",
      "Progress: 30.73%\n",
      "Progress: 30.90%\n",
      "Progress: 31.08%\n",
      "Progress: 31.25%\n",
      "Progress: 31.42%\n",
      "Progress: 31.60%\n",
      "Progress: 31.77%\n",
      "Progress: 31.94%\n",
      "Progress: 32.12%\n",
      "Progress: 32.29%\n",
      "Progress: 32.47%\n",
      "Progress: 32.64%\n",
      "Progress: 32.81%\n",
      "Progress: 32.99%\n",
      "Progress: 33.16%\n",
      "Progress: 33.33%\n",
      "Progress: 33.51%\n",
      "Progress: 33.68%\n",
      "Progress: 33.85%\n",
      "Progress: 34.03%\n",
      "Progress: 34.20%\n",
      "Progress: 34.38%\n",
      "Progress: 34.55%\n",
      "Progress: 34.72%\n",
      "Progress: 34.90%\n",
      "Progress: 35.07%\n",
      "Progress: 35.24%\n",
      "Progress: 35.42%\n",
      "Progress: 35.59%\n",
      "Progress: 35.76%\n",
      "Progress: 35.94%\n",
      "Progress: 36.11%\n",
      "Progress: 36.28%\n",
      "Progress: 36.46%\n",
      "Progress: 36.63%\n",
      "Progress: 36.81%\n",
      "Progress: 36.98%\n",
      "Progress: 37.15%\n",
      "Progress: 37.33%\n",
      "Progress: 37.50%\n",
      "Progress: 37.67%\n",
      "Progress: 37.85%\n",
      "Progress: 38.02%\n",
      "Progress: 38.19%\n",
      "Progress: 38.37%\n",
      "Progress: 38.54%\n",
      "Progress: 38.72%\n",
      "Progress: 38.89%\n",
      "Progress: 39.06%\n",
      "Progress: 39.24%\n",
      "Progress: 39.41%\n",
      "Progress: 39.58%\n",
      "Progress: 39.76%\n",
      "Progress: 39.93%\n",
      "Progress: 40.10%\n",
      "Progress: 40.28%\n",
      "Progress: 40.45%\n",
      "Progress: 40.62%\n",
      "Progress: 40.80%\n",
      "Progress: 40.97%\n",
      "Progress: 41.15%\n",
      "Progress: 41.32%\n",
      "Progress: 41.49%\n",
      "Progress: 41.67%\n",
      "Progress: 41.84%\n",
      "Progress: 42.01%\n",
      "Progress: 42.19%\n",
      "Progress: 42.36%\n",
      "Progress: 42.53%\n",
      "Progress: 42.71%\n",
      "Progress: 42.88%\n",
      "Progress: 43.06%\n",
      "Progress: 43.23%\n",
      "Progress: 43.40%\n",
      "Progress: 43.58%\n",
      "Progress: 43.75%\n",
      "Progress: 43.92%\n",
      "Progress: 44.10%\n",
      "Progress: 44.27%\n",
      "Progress: 44.44%\n",
      "Progress: 44.62%\n",
      "Progress: 44.79%\n",
      "Progress: 44.97%\n",
      "Progress: 45.14%\n",
      "Progress: 45.31%\n",
      "Progress: 45.49%\n",
      "Progress: 45.66%\n",
      "Progress: 45.83%\n",
      "Progress: 46.01%\n",
      "Progress: 46.18%\n",
      "Progress: 46.35%\n",
      "Progress: 46.53%\n",
      "Progress: 46.70%\n",
      "Progress: 46.88%\n",
      "Progress: 47.05%\n",
      "Progress: 47.22%\n",
      "Progress: 47.40%\n",
      "Progress: 47.57%\n",
      "Progress: 47.74%\n",
      "Progress: 47.92%\n",
      "Progress: 48.09%\n",
      "Progress: 48.26%\n",
      "Progress: 48.44%\n",
      "Progress: 48.61%\n",
      "Progress: 48.78%\n",
      "Progress: 48.96%\n",
      "Progress: 49.13%\n",
      "Progress: 49.31%\n",
      "Progress: 49.48%\n",
      "Progress: 49.65%\n",
      "Progress: 49.83%\n",
      "Progress: 50.00%\n",
      "Progress: 50.17%\n",
      "Progress: 50.35%\n",
      "Progress: 50.52%\n",
      "Progress: 50.69%\n",
      "Progress: 50.87%\n",
      "Progress: 51.04%\n",
      "Progress: 51.22%\n",
      "Progress: 51.39%\n",
      "Progress: 51.56%\n",
      "Progress: 51.74%\n",
      "Progress: 51.91%\n",
      "Progress: 52.08%\n",
      "Progress: 52.26%\n",
      "Progress: 52.43%\n",
      "Progress: 52.60%\n",
      "Progress: 52.78%\n",
      "Progress: 52.95%\n",
      "Progress: 53.12%\n",
      "Progress: 53.30%\n",
      "Progress: 53.47%\n",
      "Progress: 53.65%\n",
      "Progress: 53.82%\n",
      "Progress: 53.99%\n",
      "Progress: 54.17%\n",
      "Progress: 54.34%\n",
      "Progress: 54.51%\n",
      "Progress: 54.69%\n",
      "Progress: 54.86%\n",
      "Progress: 55.03%\n",
      "Progress: 55.21%\n",
      "Progress: 55.38%\n",
      "Progress: 55.56%\n",
      "Progress: 55.73%\n",
      "Progress: 55.90%\n",
      "Progress: 56.08%\n",
      "Progress: 56.25%\n",
      "Progress: 56.42%\n",
      "Progress: 56.60%\n",
      "Progress: 56.77%\n",
      "Progress: 56.94%\n",
      "Progress: 57.12%\n",
      "Progress: 57.29%\n",
      "Progress: 57.47%\n",
      "Progress: 57.64%\n",
      "Progress: 57.81%\n",
      "Progress: 57.99%\n",
      "Progress: 58.16%\n",
      "Progress: 58.33%\n",
      "Progress: 58.51%\n",
      "Progress: 58.68%\n",
      "Progress: 58.85%\n",
      "Progress: 59.03%\n",
      "Progress: 59.20%\n",
      "Progress: 59.38%\n",
      "Progress: 59.55%\n",
      "Progress: 59.72%\n",
      "Progress: 59.90%\n",
      "Progress: 60.07%\n",
      "Progress: 60.24%\n",
      "Progress: 60.42%\n",
      "Progress: 60.59%\n",
      "Progress: 60.76%\n",
      "Progress: 60.94%\n",
      "Progress: 61.11%\n",
      "Progress: 61.28%\n",
      "Progress: 61.46%\n",
      "Progress: 61.63%\n",
      "Progress: 61.81%\n",
      "Progress: 61.98%\n",
      "Progress: 62.15%\n",
      "Progress: 62.33%\n",
      "Progress: 62.50%\n",
      "Progress: 62.67%\n",
      "Progress: 62.85%\n",
      "Progress: 63.02%\n",
      "Progress: 63.19%\n",
      "Progress: 63.37%\n",
      "Progress: 63.54%\n",
      "Progress: 63.72%\n",
      "Progress: 63.89%\n",
      "Progress: 64.06%\n",
      "Progress: 64.24%\n",
      "Progress: 64.41%\n",
      "Progress: 64.58%\n",
      "Progress: 64.76%\n",
      "Progress: 64.93%\n",
      "Progress: 65.10%\n",
      "Progress: 65.28%\n",
      "Progress: 65.45%\n",
      "Progress: 65.62%\n",
      "Progress: 65.80%\n",
      "Progress: 65.97%\n",
      "Progress: 66.15%\n",
      "Progress: 66.32%\n",
      "Progress: 66.49%\n",
      "Progress: 66.67%\n",
      "Progress: 66.84%\n",
      "Progress: 67.01%\n",
      "Progress: 67.19%\n",
      "Progress: 67.36%\n",
      "Progress: 67.53%\n",
      "Progress: 67.71%\n",
      "Progress: 67.88%\n",
      "Progress: 68.06%\n",
      "Progress: 68.23%\n",
      "Progress: 68.40%\n",
      "Progress: 68.58%\n",
      "Progress: 68.75%\n",
      "Progress: 68.92%\n",
      "Progress: 69.10%\n",
      "Progress: 69.27%\n",
      "Progress: 69.44%\n",
      "Progress: 69.62%\n",
      "Progress: 69.79%\n",
      "Progress: 69.97%\n",
      "Progress: 70.14%\n",
      "Progress: 70.31%\n",
      "Progress: 70.49%\n",
      "Progress: 70.66%\n",
      "Progress: 70.83%\n",
      "Progress: 71.01%\n",
      "Progress: 71.18%\n",
      "Progress: 71.35%\n",
      "Progress: 71.53%\n",
      "Progress: 71.70%\n",
      "Progress: 71.88%\n",
      "Progress: 72.05%\n",
      "Progress: 72.22%\n",
      "Progress: 72.40%\n",
      "Progress: 72.57%\n",
      "Progress: 72.74%\n",
      "Progress: 72.92%\n",
      "Progress: 73.09%\n",
      "Progress: 73.26%\n",
      "Progress: 73.44%\n",
      "Progress: 73.61%\n",
      "Progress: 73.78%\n",
      "Progress: 73.96%\n",
      "Progress: 74.13%\n",
      "Progress: 74.31%\n",
      "Progress: 74.48%\n",
      "Progress: 74.65%\n",
      "Progress: 74.83%\n",
      "Progress: 75.00%\n",
      "Progress: 75.17%\n",
      "Progress: 75.35%\n",
      "Progress: 75.52%\n",
      "Progress: 75.69%\n",
      "Progress: 75.87%\n",
      "Progress: 76.04%\n",
      "Progress: 76.22%\n",
      "Progress: 76.39%\n",
      "Progress: 76.56%\n",
      "Progress: 76.74%\n",
      "Progress: 76.91%\n",
      "Progress: 77.08%\n",
      "Progress: 77.26%\n",
      "Progress: 77.43%\n",
      "Progress: 77.60%\n",
      "Progress: 77.78%\n",
      "Progress: 77.95%\n",
      "Progress: 78.12%\n",
      "Progress: 78.30%\n",
      "Progress: 78.47%\n",
      "Progress: 78.65%\n",
      "Progress: 78.82%\n",
      "Progress: 78.99%\n",
      "Progress: 79.17%\n",
      "Progress: 79.34%\n",
      "Progress: 79.51%\n",
      "Progress: 79.69%\n",
      "Progress: 79.86%\n",
      "Progress: 80.03%\n",
      "Progress: 80.21%\n",
      "Progress: 80.38%\n",
      "Progress: 80.56%\n",
      "Progress: 80.73%\n",
      "Progress: 80.90%\n",
      "Progress: 81.08%\n",
      "Progress: 81.25%\n",
      "Progress: 81.42%\n",
      "Progress: 81.60%\n",
      "Progress: 81.77%\n",
      "Progress: 81.94%\n",
      "Progress: 82.12%\n",
      "Progress: 82.29%\n",
      "Progress: 82.47%\n",
      "Progress: 82.64%\n",
      "Progress: 82.81%\n",
      "Progress: 82.99%\n",
      "Progress: 83.16%\n",
      "Progress: 83.33%\n",
      "Progress: 83.51%\n",
      "Progress: 83.68%\n",
      "Progress: 83.85%\n",
      "Progress: 84.03%\n",
      "Progress: 84.20%\n",
      "Progress: 84.38%\n",
      "Progress: 84.55%\n",
      "Progress: 84.72%\n",
      "Progress: 84.90%\n",
      "Progress: 85.07%\n",
      "Progress: 85.24%\n",
      "Progress: 85.42%\n",
      "Progress: 85.59%\n",
      "Progress: 85.76%\n",
      "Progress: 85.94%\n",
      "Progress: 86.11%\n",
      "Progress: 86.28%\n",
      "Progress: 86.46%\n",
      "Progress: 86.63%\n",
      "Progress: 86.81%\n",
      "Progress: 86.98%\n",
      "Progress: 87.15%\n",
      "Progress: 87.33%\n",
      "Progress: 87.50%\n",
      "Progress: 87.67%\n",
      "Progress: 87.85%\n",
      "Progress: 88.02%\n",
      "Progress: 88.19%\n",
      "Progress: 88.37%\n",
      "Progress: 88.54%\n",
      "Progress: 88.72%\n",
      "Progress: 88.89%\n",
      "Progress: 89.06%\n",
      "Progress: 89.24%\n",
      "Progress: 89.41%\n",
      "Progress: 89.58%\n",
      "Progress: 89.76%\n",
      "Progress: 89.93%\n",
      "Progress: 90.10%\n",
      "Progress: 90.28%\n",
      "Progress: 90.45%\n",
      "Progress: 90.62%\n",
      "Progress: 90.80%\n",
      "Progress: 90.97%\n",
      "Progress: 91.15%\n",
      "Progress: 91.32%\n",
      "Progress: 91.49%\n",
      "Progress: 91.67%\n",
      "Progress: 91.84%\n",
      "Progress: 92.01%\n",
      "Progress: 92.19%\n",
      "Progress: 92.36%\n",
      "Progress: 92.53%\n",
      "Progress: 92.71%\n",
      "Progress: 92.88%\n",
      "Progress: 93.06%\n",
      "Progress: 93.23%\n",
      "Progress: 93.40%\n",
      "Progress: 93.58%\n",
      "Progress: 93.75%\n",
      "Progress: 93.92%\n",
      "Progress: 94.10%\n",
      "Progress: 94.27%\n",
      "Progress: 94.44%\n",
      "Progress: 94.62%\n",
      "Progress: 94.79%\n",
      "Progress: 94.97%\n",
      "Progress: 95.14%\n",
      "Progress: 95.31%\n",
      "Progress: 95.49%\n",
      "Progress: 95.66%\n",
      "Progress: 95.83%\n",
      "Progress: 96.01%\n",
      "Progress: 96.18%\n",
      "Progress: 96.35%\n",
      "Progress: 96.53%\n",
      "Progress: 96.70%\n",
      "Progress: 96.88%\n",
      "Progress: 97.05%\n",
      "Progress: 97.22%\n",
      "Progress: 97.40%\n",
      "Progress: 97.57%\n",
      "Progress: 97.74%\n",
      "Progress: 97.92%\n",
      "Progress: 98.09%\n",
      "Progress: 98.26%\n",
      "Progress: 98.44%\n",
      "Progress: 98.61%\n",
      "Progress: 98.78%\n",
      "Progress: 98.96%\n",
      "Progress: 99.13%\n",
      "Progress: 99.31%\n",
      "Progress: 99.48%\n",
      "Progress: 99.65%\n",
      "Progress: 99.83%\n",
      "Progress: 100.00%\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "'done'"
      ]
     },
     "execution_count": 40,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "dump_discrete_log(L, P)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "759d4819-89b1-4b77-897b-a0bc175564ef",
   "metadata": {},
   "outputs": [],
   "source": [
    "def load_csv_as_matrix(filename):\n",
    "    \"\"\"\n",
    "    load a .csv file which is a matrix of integers\n",
    "    \"\"\"\n",
    "    import csv\n",
    "    with open(filename, newline='') as csvfile:\n",
    "        reader = csv.reader(csvfile)\n",
    "        matrix = [list(map(int, row)) for row in reader]\n",
    "    return Matrix(matrix)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "0ce9ef07-5d27-47a9-9291-c50b8e0db817",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "1723872719529775255168233477739808705550682172067583351531515002896576794274699331456573828897879671480606241901141153738537923622099791332217469260702699225221018982963916207228\n"
     ]
    }
   ],
   "source": [
    "filename = f\"data/discrete_log_eigenvector_matrix/discrete_log_dft_symmetric_group_finite_field_eigenvector_matrix_n={n}_q={q}.csv\"\n",
    "dlog_P = load_csv_as_matrix(filename); print(dlog_P[1,0])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 115,
   "id": "f8659c2a-722c-476d-8073-8c25f2628be8",
   "metadata": {},
   "outputs": [],
   "source": [
    "#compute log1p of dlog_P, which makes the magnitudes much smaller\n",
    "dlop_P_numeric = np.array(dlog_P, dtype=np.float64)\n",
    "log1p_dlog_P = np.log1p(dlop_P_numeric)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 125,
   "id": "00a2ebf6-1564-4ee3-903b-01015ff42ded",
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "from matplotlib.colors import ListedColormap, BoundaryNorm\n",
    "\n",
    "def plot_log1p_dlog(F, M, path, title):\n",
    "    \n",
    "    M_numeric = np.array(M, dtype=np.float64)\n",
    "    log1p_M = np.log1p(M_numeric) #apply log1p transformation\n",
    "    pos_log1p_M = log1p_M[(log1p_M != -np.inf) & (log1p_M != 0)] #mask out the -np.inf values and zero values for normalization\n",
    "\n",
    "    cmap = plt.cm.hsv  # Use the HSV colormap\n",
    "    num_colors = min(F.order(), 256)  # Limit the number of colors\n",
    "    new_colors = np.vstack(([0, 0, 0, 1], cmap(np.linspace(0, 1, num_colors))))\n",
    "    custom_cmap = ListedColormap(new_colors)  # Create a new custom colormap\n",
    "\n",
    "    # Normalize: exclude -np.inf values for the min/max scaling\n",
    "    norm = BoundaryNorm([np.min(pos_log1p_M), *np.linspace(np.min(pos_log1p_M), np.max(pos_log1p_M), num_colors)], custom_cmap.N)\n",
    "\n",
    "    # Plot the matrix with the custom colormap\n",
    "    plt.imshow(log1p_M, cmap=custom_cmap, norm=norm, interpolation=\"nearest\")\n",
    "    plt.title(title, fontsize=16)\n",
    "    plt.colorbar()\n",
    "    plot_title = path + title.replace(' ', '_') + '.png'\n",
    "    plt.savefig(plot_title, dpi=300, bbox_inches=\"tight\")\n",
    "    plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 126,
   "id": "af605ddb-d4ab-4f68-9fb9-d5fec81f5255",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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VekL58MPPyyKi4ul+z7++GPpudm6davRfpqTZN7qTiS2eXl50nvz+fPnpfWNndjq1fVecPbsWSnxnT17tigvLxdC1Lz+Xn/9dem19OOPPxrE6R9XtVot7r//foNz49q1a2Yfm4ODgwgICBDff/+9dN/Ro0elLzvLly+vFfvggw8afMExdzHm2rVrokOHDqJly5bSa05JYnvgwAGhVquFSqUS6enpUsJWUlIixo8fL73fGUsALXlcxowYMUIANYM9d+q4SktLpS/GTz/9tCgqKpLuO3TokOjatasAINLT0w3i9AMerVq1Env37pXWnzlzRvTs2VNq79bXrf4zuqGLuYNoeg888IAAIP7yl7/Uu605iW15ebk0CLFv3z6j2+zdu1f6stmQLyYWSWyLioqkD/JNmzbVivnPf/4jAAiVSiWOHz8urV+xYoUAIDp16iQqKytrxenfSG5HYtujRw+j9+sTxr/97W/1tilETZKpH60y9oH322+/SaOe33zzjdE+NzSx/frrrwVQM1r2xx9/1Lr/6aeflh43cxJbucewbds26dvwlStXasXpP0AbktjKPZdMHZs5tmzZIgCIiIgIo+fhuXPnhJubm3BzczMYLTWVfL3wwgsCgBg6dKjR9vT9vDWBkNsP/ePs7Ows8vPza8WNHTtWADV/AbnZY489JgAYHW369NNPBQAxYMAAg/Vz5swRAMQzzzxj9Nj+9a9/CQBi8ODB0rqysjLprwtr1qwxGner+hJbueeJPmEEIN5++22j+54xY4YAYJD0yrVy5UoBQPj5+RlNfuLj440+zk09sb1y5Yr0ON6cAN6cmJpa6nq/uxOJ7fz58wVQM6JvzH333ScAiClTphis1z+uWq1WlJaWNrhPNx/bxo0ba93/zjvvCADi/vvvr3Wf/nOroYsx+uNfsWJFrb7JSWwfeeQRAUA89NBDte67du2a9JdJY4mtJY/rVgUFBdIXx+++++6OHZf+eRwzZozR/f74449CpVIZjEJWVVVJgwLGRi+PHTsmHcutr9ub39MasjRkBH379u1SXHZ2dr3bm5PYFhYWSvu8edDlZocPH5a2MfYXf1MsUmO7Z88elJWVISgoCKNHj651f0REBPr06QMhBLKysqT1+v9PmTLFaI3P7ZxzLj4+vs7127dvN2s/R44cQX5+PpycnBAXF1fr/rvuukuaFicz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L774Yr3xHTt2xI8//ojDhw8bvf/YsWOoqKiAWq02+SfqhtJ/SJr64nv8+HFF+6+urkZGRgYAGH3fbuo6duwInU6Hw4cPG/0w1Z83fn5+DRolup10Oh1Onz5tkX3VNSd5VVWVlNib+vPwzVq0aAFfX18UFhbiv//9r9EEUP+Zaowlj0tPCIG1a9cCkD93rZLj6tKlC/bt24eff/7Z7MT25s/0yspKo0mcqc+N69evm/UXtFuZmrv4xo0bePjhh7Fjxw60a9cOWVlZJsu2lIqMjMSpU6ewd+9eo7+/0B9XZGRkg/ZrkYK3/v37w8XFBfn5+UavEPHDDz/gu+++g0qlwrBhw6T1+v9/8sknqKqqqhVnaoJxSzA1ybp+vX4y4vp07twZQUFBuH79Oj744INa9587d0661O/w4cNl9taQvm+fffaZ0ZHZhk4gL/cY9M/fF198YfTbpJznT+65pIT+hxzvv/++QUG7XPp+mTp+U+st3Q9z6UdmP/30U2zYsAFVVVUYNWqU0S84+j6mpaWZvf9+/frB29sbZ8+eNfpDNWP09VSmSj8a4zxRq9XSD3zOnTtnVoz+9bJy5Uqjz+k777wDoOYxUvonbj19HbBOp6t132+//Wb2D2NNWbBgAX788Uc4ODjg2WefVbSvxqB/TtLT043er39OLPV+bQmW+JFVWlqaye2+/fZbADX1v/p15l5QQf/60l+U52bl5eX46KOPbutx3So7OxtnzpxRfAlducelf4985513zC4r6dy5M+666y6UlpYavcjNr7/+avLiUZb88ZgQAtOmTcOWLVsQEBCAr7/+GgEBAWYdgxz6x+qf//xnrRyipKREujBDg0ek6540obb65rG96667DCYRP378uDTh/Pjx4w1ibp7H9oknnpCmw6murhbLli277fPYPvfccwbz2L744osCqJnL9ea5QOujnwPWw8NDfP3119L6goICae7R3r17m+yznHlse/bsKYDa88du2LBBODo6yp7HtiHHUFVVJbp37y6AmrntLl++LN23adMm4eTkpGge24acS3UdW32qqqpE7969BVAzRdexY8cM7r9+/br497//LU1mrmfOPLYLFy40mMd2/vz5dc5jK6cfpqZvqq+fekeOHBFAzYUk9BcG2Lx5s9Ft8/PzRatWrQRQM+XXzc+5EDVz0n744Ye1XodLly4VQM0FCtatW2cwOfmlS5dEamqqwTy2Op1OABDBwcHi6tWrRvsi5zwxZxqtl156SXzwwQe1ji0vL08EBAQIAOKjjz4yGX+zm+exHT9+vCgpKZHu+8c//iG9Ti05j61+ihx7e3vx1VdfSevPnTsnBg4caHI+zLqm+6qqqhL79u0T999/v/Se/MEHH9TaztrmsU1MTDSYx1Y/IX9d89g29P1a7+YLNBij5DlXSsl0Xz/88IOws7MTKpVKLF++XHptl5aWikmTJtU5j+3tMH36dAFAjBs3TtF+5B5XSUmJaNeunQAgJk6cWGtO/JKSErFhw4ZaUybq57H19vY2mJ4sPz9fhIWFmXzdWpL+IhDe3t7i8OHDsvdj7gUabty4Ic05Pm7cOOm9vrS0VIwbN04ANRdZuXnec3NY9MpjgwYNkt6UunTpInr06CHNvdajRw+TVx7Tz5Hq6ekpIiIipA+Pt99+WwA1k/Peas+ePcLLy0ta9Emwi4uLwfqbJ5e/uf/6K495eXmJiIgI4e3tLbX16aefNugxufWqXR06dDC4aldQUJDRq3YpeaP8+eefpQTD2dlZuooTYJkrj5l7DDdfeczFxUWEh4dLVxOZNWuW1I9bn4e6yD2X5Ca2QtR86Ou/LOiPPzIyUnTp0kV6DPz8/Axi6koYP/74Y+kqMz4+PtI5ZmdnJ5YsWSIA41dTkdMPpYmtEELcc889UpstWrSQPuyN2bNnj/R6cXBwEN26dRORkZGiXbt20jHfmlBWV1eLv/71r1Ib+quxtW3bVnpebz5Xq6qqpKuPeXl5iT59+oioqCgxe/ZsaRs554k5ycPo0aOl94IOHTqIe++9V3To0EFqZ9CgQUbnGTZl//790sUwXF1dRXh4uMGV2l5++eVaMUqTHP2Hu/7Lgf4KRp06dRKzZ8+uM7G9+Wpnffv2FaGhocLd3V3an4+Pj/jiiy+MtmsNia0QNVce07+eWrZsKSIiIqR5Se3s7MT7779fK4aJrWnJycnS8xYQECDCw8OFu7v7Hb/y2LVr16TXmqkv5w0h97iOHDkigoODpfOpc+fOIjIyUnTs2FF6f7r1gkeVlZUiJiZGaq9Tp06iZ8+eta489uqrryo+LmP27dsntR0YGNjgC2TcnHe5ublJnw83rzeWW+Xl5YmWLVtKOWBYWJj0HLZq1UocOnSowcdiscRWCCEqKirE0qVLRXh4uHB1dZUut7lo0SKToy5CCJGTkyNGjhwpPD09haurq4iIiBCffvqpKC0tlQ72VuZMBG6snzf3f9u2bWLgwIHC3d1duLm5icGDB5s1+bAx1dXVYs2aNWLAgAHCw8NDaDQaERISIp577jmTEwsrfaP89ddfxaRJk4SXl5dwcnIS3bp1E++++66orq5ucGIr9xiEqBkhmzBhgkE/9FdV0SdAt45+1UfOuaQksRWiZkR02bJlYuDAgaJly5bC0dFRBAYGiv79+4tXXnml1jfY+hLGrKwsER0dLdzc3IS7u7uIiooSmZmZ4ueff5YSL0v0wxKJrT7ZBiCmT59ucju9wsJCMX/+fNGjRw/h5uYmnJ2dRYcOHURsbKxYtmyZwVXTbvbVV1+Jv/zlL8LHx0c4OjqKu+66SwwePFgsW7asVjL9yy+/iAcffFD4+vpKHwa3fog09DwxJ3nQ6XTihRdeEJGRkUKr1Ur9jIqKEmvWrGlQUqt36tQpMWPGDNGmTRvh6OgoWrZsKWJiYgxGVBvaz7pUVlaKV199VbRv317qf0JCgrh8+XK9F2i4eXFwcBDe3t4iNDRUPProo+If//hHnZd8tZbEVoiaD9TJkycLf39/4eDgIHx8fMSYMWNMXgWJiW3dPv/8cxEZGSmcnZ1Fy5YtxX333Sd0Ol29x21JGzZskL4My7mErjFyj6u4uFi88cYbIjIyUvo8bdu2rRg8eLB46623jOZQ5eXlYtGiRSIkJEQ4OjoKf39/MX36dPH777+LuXPnCsD4RZQswdycytRr05w4U6+dM2fOiCeeeELcdddd0vtVXFyc0QsOmUP1vw41STk5OQgPD0ePHj2Qm5trkX22bdsWp0+fxsmTJ836VTPJd/HiRXh7e6NFixYmp/NojjZu3IgHH3wQo0ePxubNmxu7O0REt83OnTsxaNAgREVF3fF5WG3JqFGj8O9//xubNm2q9cNeMmSZ2dJvk1WrVgGo+WEFWR/981fflDHNDc9rIiIy12+//YasrCyo1WpZV1Frbho9sf3222+xfv16gytMVFZWIjU1FcuXL4ednZ3FLv1HlpeXl4eMjAyD+feEEPjkk0+wYMECAMBTTz3VWN1rNBs3bsTWrVsNZvsoKyvDvHnz8NVXX8HV1VX2VDRERGR7Fi1aJF2FTO/o0aMYPXq0NLessanHyJBF5rFV4vTp03jsscfg4OCA4OBgeHh44JdffpHmP01JSTF7vky68y5evIgZM2YgPj4ebdq0gZeXF3799VdpGrIZM2Zg1KhRjdzLOy8vLw+vvPIKnJyc0L59e2g0Ghw5cgTXrl2DWq3G+++/zzcoIiKSfPDBB1iwYAG8vb3Rtm1bFBUVSYluu3btpKnoqG6NPmI7YMAAzJw5Ex07dsQff/yB3NxcODk5YdSoUdi+fTteeOGFxu4i1aFLly6YN28eunXrhqKiIhw8eBBCCAwZMgTr1683OgdgczB69GhMnz4dgYGByM/PR15eHlq2bInx48fju+++MzoZNRERNV8LFizA8OHDodFo8PPPP+Ps2bPo2rUr5s+fjx9++AF33XVXY3fRKjTpH48REREREZmr0UdsiYiIiIgsodFrbImqq6tx7tw5uLu7Q6VSNXZ3iIiogYQQKCkpQUBAAOzsbt+Y2fXr11FRUaF4P46OjnBycrJAj6ipYWJLje7cuXMIDAxs7G4QEZFC+fn5aN269W3Z9/Xr1xHs7IoCVCvel4eHB/z9/WFnZ4eEhAQkJCRYoIfUFDCxpUbn7u7e2F0gIiILuJ3v5xUVFShANfLRFh4KKimLUY3A4lOcmcZGMbGlRsfyAyIi23An3s89YKcosdXT6XTw8PCwQI+oKeGPx4iIiMiKOABwVLA4AAAiIiLQpUsXvPfee3e4/3Q7MbEli1m2bBmCg4Ph5OSEsLAw7N69u7G7RERENkdJUqtfyFYxsSWL2LBhAxITEzF//nwcPHgQAwYMQGxsLM6cOdPYXSMiIpviYIGFbBUv0EAWERkZiV69emH58uXSus6dO+OBBx5ASkqKwbbl5eUoLy+XbhcXF3NWBCIiG1BUVHTb6laLi4vh6emJIvSCB9Ty94MqeOLAbe0rNR6O2JJiFRUVyMnJQUxMjMH6mJgY7Nu3r9b2KSkp8PT0lBYmtUREZD7W2JJpTGxJsQsXLqCqqgp+fn4G6/38/FBQUFBr+xdffBFFRUXSkp+ff6e6SkREVo81tmQap/sii7l1mhchhNGpXzQaDTQazZ3qFhERUS2c7ss2ccSWFPP29oZara41OltYWFhrFJeIiEgZy4zYshTBNjGxJcUcHR0RFhaGrKwsg/VZWVno27dvI/WKiIhsk2VqbMk2sRSBLGLOnDmYMmUKwsPD0adPH2RkZODMmTN46qmnGrtrREREtbAUwTYxsSWLGD9+PC5evIhXX30V58+fR2hoKLZu3Yo2bdo0dteIiMimKB11rfntR0REBNRqNRISEpCQkGCRnlHjY2JLFhMfH4/4+HjZ8UVIggecZES2lN1mjauyI5fgGUUtr1cQe/Af8mP/PkVBwwAqFMTOQxdljeNL2ZHbEaKo5REfyo8V0zvIjlUJBU82AKEqkh37LUYoansQ7pEdq3ozV37Dc8fLjwXwqWqD7NgJeE1R25exQHZsS/xdUduqKzLf04qLgSBPRW2bzzKJLdkmJrZERETU7LAUwTbxx2NERERkRXiBBjKNI7ZERERkRXiRBTKNiS0RERE1OyxFsE0sRSAiIiIrwgs0kGkcsSUiIiIroq+xlUtYqiPUBDGxJSIiIiviAGXTfVVbqiPUBDGxJSIiomaHNba2iTW2REREZEVYY0umccSWiIiIrIjSGluWItgyJrZERETU7LAUwTaxFIGIiIisCEsRyDSO2BIREZEVUXrlMZYi2DKO2BIREVGzo9PpcPjwYSQkJDQoLiUlBSqVComJidI6IQSSkpIQEBAAZ2dnREdH49ChQwZxBQUFmDJlCrRaLVxdXdGrVy98/vnn9ba3bNkyBAcHw8nJCWFhYdi9e7fB/bezbWvExJaIiIisiP7HY3KXmjlw5ZQi6HQ6ZGRkoHv37gbrlyxZgtTUVKSnp0On00Gr1WLYsGEoKSmRtpkyZQqOHj2KLVu2IC8vD2PHjsX48eNx8OBBk+1t2LABiYmJmD9/Pg4ePIgBAwYgNjYWZ86cue1tWyuVEIKX4KBGVVxcDE9PTxQhEh4yqmNU4nFF7S9WyY+fh2OK2lZmtILYqQrb3q4gNkJh260UxAYqank7JsmOfVZBuy0VxALAywpitQrb7ooy2bEZcJEdG6+4511kR87BDkUtp+IeBdFuitpWLd8jL/BaMTDHE0VFRbftB1l/flYshQec5e8H1+CJ2ejYsSPUajUSEhLMGrUtLS1Fr169sGzZMixatAj33HMP0tLSIIRAQEAAEhMT8fzzzwMAysvL4efnh8WLF2PGjBkAADc3NyxfvhxTpkyR9unl5YUlS5Zg+vTpRtuMjIxEr169sHz5cmld586d8cADDyAlJeW2tm2tOGJLREREzc4333yD/fv3Y8qUKSguLkZ5eXmd2yckJGDkyJEYOnSowfqTJ0+ioKAAMTEx0jqNRoOoqCjs27dPWte/f39s2LABly5dQnV1NdavX4/y8nJER0cbba+iogI5OTkG+wWAmJgYab+3q21rxh+PERERkRVR+uOxKgBAYKDhX3AWLlyIpKQkoxHr16/HgQMHoNPpat1XUFAAAPDz8zNY7+fnh9OnT0u3N2zYgPHjx8PLywv29vZwcXHBpk2b0L59e6NtXrhwAVVVVUb3q2/zdrVtzZjYEhERkRVReoGGGwCADh06QK1WIy4uDnFxcdBoNEa3zs/Px+zZs5GZmQknJyeTe1WpVAa3hRAG615++WVcvnwZX3/9Nby9vbF582Y89NBD2L17N7p16yZ7v7ezbWvExJaIiIisiNIR25rE1s7ODnZ2dnBycqqzLjgnJweFhYUICwuT1lVVVWHXrl1IT0/H0aNHAdSMnvr7+0vbFBYWSiOpJ06cQHp6On7++Wd07doVANCjRw/s3r0b7733HlasWFGrXW9vb6jVamlU1th+tVrtbWnbmrHGloiIiJodc6f7GjJkCPLy8pCbmyst4eHhmDx5MnJzc9GuXTtotVpkZWVJMRUVFcjOzkbfvn0BAGVlNT+itLMzTLvUajWqq43Pq+vo6IiwsDCD/QJAVlaWtN/g4ODb0rY144gtERERWREH6Kfskh9fM92XObMiuLu7IzQ01GCdq6srvLy8pPWJiYlITk5GSEgIQkJCkJycDBcXF0yaVDOTSqdOndChQwfMmDEDb731Fry8vLB582ZkZWXh3//+t7TfIUOGYMyYMZg5cyYAYM6cOZgyZQrCw8PRp08fZGRk4MyZM3jqqacAQJpP1xJt2womtkRERGRFlNbYVlqqI5J58+bh2rVriI+Px+XLlxEZGYnMzEy4u7sDABwcHLB161a88MILGDVqFEpLS9GhQwesXr0a9913n7SfEydO4MKFC9Lt8ePH4+LFi3j11Vdx/vx5hIaGYuvWrWjTpo3F27YVnMeWGh3nsZWL89g2HOexbSjOY9twnMf2ds9j+zk84Cp/P7gKTzx4W/tKjYc1tkRERGRFlFx17M8fnsm58hg1fSxFICIiIiuidFYEy5ciUNPBxJaIiIiaHZ1Ox1IEG8RSBCIiIrIi+h+PyV3+nBWBpQi2hyO2REREZEWUliIoiaWmjoktERERNTssRbBNTGypCRkNwPR1uE0Rqu8UtttRQewchW3Hy478HYdlx/5ddmSNN/Ctguh/KWw9UXZkWwXTdQHAKTwjO7afgkfdTdH0T8AS5MqO/UBRy8Bn+Kvs2ARxSnZsgaqg/o3qoGRCvFSMU9Q28IbsyDUYoahl8dd1suKKUQZPRS03hGVGbM29QANZFya2REREZEWUXqBByVXLqKljYktERETNDksRbBNnRSAiIiIrwgs0kGkcsSUiIiIrwlkRyDQmtkRERGRFWGNLpjGxJSIiomaHNba2iTW2REREZDUE7BUvAGtsbRVHbImIiMhqVMMe1QrSFyWx1PTx2SUiIqJmh6UItomlCERERGQ19CO2ShaApQi2iiO2REREZDVYikB14bNLREREzQ5LEWwTSxGIiIjIatyAneIFYCmCreKILREREVmNiv8tSuLJdjGxJSIiomaHpQi2iYktNRkViEMFGv4mY4cbitq1R3cF0RMVtQ1EyI50VdDqG+igIBpQdtydFbX8LzwpO3aVopYBIFF2pLungmZfzlUQDIjnDimIXqKobaBUdqRQ9VLQrtJz/HPZkQ9ho6KW31AQ+ygGK2obKJAZV6WwXfNVQtmoa+X//o2IiIBarUZCQgISEhIs0DNqCpjYEhERkdVgKQLVhYktERERWQ0mtlQXJrZERETU7LDG1jZxui8iIiKyGvoaW7nLzTW2cqb7SklJgUqlQmJiorROCIGkpCQEBATA2dkZ0dHROHToz9r2U6dOQaVSGV0+++yzOtuKiIiAu7s7fH198cADD+Do0aO1tjty5Ajuv/9+eHp6wt3dHb1798aZM2dqbSeEQGxsLFQqFTZv3tyg47YWTGyJiIjIaihJapWWMeh0OmRkZKB7d8MfHS9ZsgSpqalIT0+HTqeDVqvFsGHDUFJSAgAIDAzE+fPnDZZXXnkFrq6uiI2NNdlednY2EhISsH//fmRlZeHGjRuIiYnB1atXpW1OnDiB/v37o1OnTti5cyd+/PFHLFiwAE5OTrX2l5aWBpVKpeARaPpYikBERETNTkNLEUpLSzF58mSsXLkSixYtktYLIZCWlob58+dj7NixAIDVq1fDz88P69atw4wZM6BWq6HVag32t2nTJowfPx5ubm4m29y2bZvB7VWrVsHX1xc5OTkYOHAgAGD+/Pm47777sGTJnzOYtGvXrta+fvzxR6SmpkKn08Hf39/s47Y2HLElIiIiq1FpgQUAwsLC0KlTJ7z99tsoLi5GeXl5ne0mJCRg5MiRGDp0qMH6kydPoqCgADExMdI6jUaDqKgo7Nu3z+i+cnJykJubi+nTpzfo2IuKigAArVq1AgBUV1fjq6++QseOHTF8+HD4+voiMjKyVplBWVkZJk6ciPT09FoJtq1hYktERERWw1I1tsePH8fRo0cxd+5ceHp6IiUlxWSb69evx4EDB4xuU1BQM/evn5+fwXo/Pz/pvlt9+OGH6Ny5M/r27WveQaNmZHjOnDno378/QkNDAQCFhYUoLS3FG2+8gREjRiAzMxNjxozB2LFjkZ2dLcU+88wz6Nu3L0aPHm12e9aKpQhERETU7OTn5xuUImg0GpPbzZ49G5mZmUbrVvVurV0VQhitZ7127RrWrVuHBQsWNKi/M2fOxE8//YQ9e/ZI66qrqwEAo0ePxjPPPAMAuOeee7Bv3z6sWLECUVFR2LJlC3bs2IGDBw82qD1rxRFbIiIishqW+vHYkCFD0Lt3b/zjH/+Ah4eHycQ2JycHhYWFCAsLg729Pezt7ZGdnY133nkH9vb20kjtraOzhYWFtUZxAeDzzz9HWVkZHn30UbOPedasWdiyZQu+/fZbtG7dWlrv7e0Ne3t7dOnSxWD7zp07S7Mi7NixAydOnECLFi2k/gPAuHHjEB0dbXYfrAVHbImIiMhqVABwUBjfEEOGDEFeXp7BusceewydOnXC888/j3bt2kGr1SIrKws9e/asaaOiAtnZ2Vi8eHGt/X344Ye4//774ePjU2/bQgjMmjULmzZtws6dOxEcHGxwv6OjIyIiImpNAfbLL7+gTZs2AIAXXngBTzzxhMH93bp1w9///neMGjWq/gfAyjCxJSIiombH3FkR3N3dpZpWPVdXV3h5eUnrExMTkZycjJCQEISEhCA5ORkuLi6YNGmSQdzx48exa9cubN261WhbQ4YMwZgxYzBz5kwANT9YW7duHb788ku4u7tLo8Kenp5wdnYGADz33HMYP348Bg4ciEGDBmHbtm3417/+hZ07dwIAtFqt0R+MBQUF1UqUbQETWyIiIrIa+h+PKYkHai7QoFarkZCQgISEBEV9mjdvHq5du4b4+HhcvnwZkZGRyMzMhLu7u8F2H330Ee666y6DGRRuduLECVy4cEG6vXz5cgCoVTKwatUqTJs2DQAwZswYrFixAikpKXj66adx9913Y+PGjejfv7+iY7JWTGyJiIjIalRAWfKiJCnW04+G6qlUKiQlJSEpKanOuOTkZCQnJ5u8/9SpUwa3hRBm9efxxx/H448/bta2DdmvNWJiS02GI/rBEeoGx5XisKJ23T2Py44VRa8pahtYJjtyLxp2Gcib/Qj5xwwA+UL+lWvSxyt8Qz0iP1TkDVbWNlJlR1YW1V9PZ4r9cw2b6/JWa/Cp7NjHnVYrartKftMQY+5X0PJEBbEA8B/ZkZ+h9g+GGuIhyJ88/zO8r6htwPTFAupWAuBuhW2bpykkttR0MbElIiKiZqehVx4j68DpvoiIiMhqWOoCDREREejSpQvee0/+X7+o6eGILREREVmNCkBG0ZphPNkuJrZERETU7LAUwTaxFIEUS0pKgkqlMliMzZlHRESklKWuPMZSBNvEEVuyiK5du+Lrr7+WbqvVSv5QREREZNwN/FknKzeebBcTW7IIe3t7s0dpy8vLUV5eLt0uLi6+Xd0iIiIyiqUItomlCGQRx44dQ0BAAIKDgzFhwgT8+uuvJrdNSUmBp6entAQGBt7BnhIRkTVjKQLVhYktKRYZGYk1a9Zg+/btWLlyJQoKCtC3b19cvHjR6PYvvvgiioqKpCU/P/8O95iIiKyVpRJbsk0sRSDFYmNjpf9369YNffr0Qfv27bF69WrMmTOn1vYajQYajeZOdpGIiMgASxFsE0dsyeJcXV3RrVs3HDt2rLG7QkRENoYXaKC6MLEliysvL8eRI0fg7y//eudERETGsBSB6sJSBFJs7ty5GDVqFIKCglBYWIhFixahuLgYU6dObeyuERERGcVSBNvEEVtS7LfffsPEiRNx9913Y+zYsXB0dMT+/fvRpk2bxu4aERHZGM6KQHXhiC0ptn79eovsZzB0UKPh356/h7KrnIki+fGlyFXUthu2yY4djgwFsV/IjgUAFQbKji3/p7I/BDriNdmxi7BDUdsLPOXHiiIlf8G4V0Es8Cjulx279/oWRW2/P2ai7NgTkN92eyTIjq2h5DxNU9TyZ9irKF6JUOyRFVcFVwv3xLRKACqF8WS7mNgSERGR1VBaI8saW9vGxJaIiIiaHdbY2ibW2BIREZHVYI0t1YUjtkRERGQ1lNbIssbWtjGxJSIiomaHpQi2iaUIREREZDVYikB14YgtERERWQ2WIlBdmNgSERFRs8NSBNvEUgQiIiKyGixFoLpwxJaIiIisRiUAoSD+hqU6Qk0SE1siIiJqdliKYJtYikBERERWo7FLEVJSUqBSqZCYmCitE0IgKSkJAQEBcHZ2RnR0NA4dOlQr9rvvvsPgwYPh6uqKFi1aIDo6GteuXTPZVlJSElQqlcGi1Wql+ysrK/H888+jW7ducHV1RUBAAB599FGcO3fOYD8ZGRmIjo6Gh4cHVCoVrly50qBjtiZMbImIiMhqWCqxlUOn0yEjIwPdu3c3WL9kyRKkpqYiPT0dOp0OWq0Ww4YNQ0lJibTNd999hxEjRiAmJgb/+c9/oNPpMHPmTNjZ1Z2Kde3aFefPn5eWvLw86b6ysjIcOHAACxYswIEDB/DFF1/gl19+wf3332+wj7KyMowYMQIvvfSSgqO3DixFICIiIqtRCaBaQXyVzLjS0lJMnjwZK1euxKJFi6T1QgikpaVh/vz5GDt2LABg9erV8PPzw7p16zBjxgwAwDPPPIOnn34aL7zwghQbEhJSb7v29vYGo7Q38/T0RFZWlsG6d999F/feey/OnDmDoKAgAJBGl3fu3Gn28VorjtgSERFRs/PNN99g//79mDJlCoqLi1FeXl7n9gkJCRg5ciSGDh1qsP7kyZMoKChATEyMtE6j0SAqKgr79u0DABQWFuL777+Hr68v+vbtCz8/P0RFRWHPnj319vPYsWMICAhAcHAwJkyYgF9//bXO7YuKiqBSqdCiRYt6922LOGJLTcavkPtNy01RuyqfXNmx4o8Fytr+8DXZsVumy2/XEU/KDwZQpSqWHWuHwYraVuJl3KssvmiiguhWsiP/hmcUtAu85iM/XvyxQ1Hb1XhfduwcBe1+iR4KogHVV/JfI2Jk/clK3drLjnxB4Wt7rcy4UgD9FbVsvgoAagXx+hHbwMBAg/ULFy5EUlKS0Zj169fjwIED0Ol0te4rKCgAAPj5+Rms9/Pzw+nTpwFASkaTkpLw1ltv4Z577sGaNWswZMgQ/PzzzyZHbiMjI7FmzRp07NgRv//+OxYtWoS+ffvi0KFD8PLyqrX99evX8cILL2DSpEnN9odxTGyJiIjIalgqse3QoQPUajXi4uIQFxcHjUZjdPv8/HzMnj0bmZmZcHJyMrlflUplcFsIIa2rrq4pnpgxYwYee+wxAEDPnj3xzTff4KOPPkJKSorRfcbGxkr/79atG/r06YP27dtj9erVmDPH8KtfZWUlJkyYgOrqaixbtsz0A2DjmNgSERFRs5OTk2PWqGZOTg4KCwsRFhYmrauqqsKuXbuQnp6Oo0ePAqgZufX395e2KSwslEZx9eu7dOlisO/OnTvjzJkzZvfZ1dUV3bp1w7FjxwzWV1ZW4uGHH8bJkyexY8eOZjtaC7DGloiIiKxIJZTNiFD5v/2YO93XkCFDkJeXh9zcXGkJDw/H5MmTkZubi3bt2kGr1Rr8iKuiogLZ2dno27cvAKBt27YICAiQkmC9X375BW3atDH72MvLy3HkyBGDBFqf1B47dgxff/210RKF5oQjtkRERGQ1KqBsVK6hMyq4u7sjNDTUYJ2rqyu8vLyk9YmJiUhOTkZISAhCQkKQnJwMFxcXTJo0CUBNmcJzzz2HhQsXokePHrjnnnuwevVq/Pe//8Xnn38u7XfIkCEYM2YMZs6cCQCYO3cuRo0ahaCgIBQWFmLRokUoLi7G1KlTAQA3btzAgw8+iAMHDuDf//43qqqqpJrfVq1awdHREUDNaHJBQQGOHz8OAMjLy4O7uzuCgoLQqpX8+v+miIktERERNTuWvPLYvHnzcO3aNcTHx+Py5cuIjIxEZmYm3N3dpW0SExNx/fp1PPPMM7h06RJ69OiBrKwstG//548FT5w4gQsXLki3f/vtN0ycOBEXLlyAj48Pevfujf3790ujvL/99hu2bNkCALjnnnsM+vTtt98iOjoaALBixQq88sor0n0DBw4EAKxatQrTpk2zyGPQVKiEEEouuUykWHFxMTw9PdESRbBDw99kLqCDovZVPsdlx1rvrAjKDIOSWRHuU9i6EkqmZgcAJbMimF9Hd6u/4e8K2gVe85Efq3xWhHDZsWNkvB/ofYkM2bFAY8+KcFh2pNJZEeSe4aUoRn94oqio6LbVdyr9rNCrRjEuwxMdO3aEWq1GQkICEhISLNhTakwcsSUiIiKrcQOAqt6tTONonm1jYktERETNjiVLEajp4KwIREREZDWUzIigXwDzZ0Ug68IRWyIiIrIa1WpApaAWQQj8eZUGsjlMbImIiKjZYSmCbWIpAhEREVmNanvlC8BSBFvFEVsiIiKyGtX2FihFKLdYd6iJYWJLREREVkPYA4LzfZEJTGypyfgZkDnldoSidsUfU2XHXsZqRW1rFFxkob+CdlviXgXRAHBJQewbilpW/fCg7FgRPk5Z2yOfkd/2V/IvsvAqBsuOBQBHRRdZUDbhvx36yY79UtGlRP6lIBYQI39XEC3/Ags1BsmOfAPKzvEl2Cgr7rqiVhsHa2xtE2tsiYiIyHo4oObyiXIXh5rdsMbWNnHEloiIiKyHI5QNy1VbqiPUFDGxJSIiomaHpQi2iaUIREREZD2UlCHoF7AUwVZxxJaIiIishwMAtYJ4XnXMpjGxJSIiomaHpQi2iaUIREREZD0cLLCApQi2iiO2REREZD0cwVIEMomJLRERETU7LEWwTSxFICIiIuvBCzRQHThiS0RERNbDEcqyFw7p2TQmtkRERGQ9mNhSHZjYEhERUbPDGlvbxO8tREREZD1YY0t14IgtERERWY+bklNZVJbqCDVFTGypyXDDIrhBIyPyiKJ2b2Cf7NjDiloGrqOD7NgvcFx27FgMlx0LAKrlr8mOXf1XRU1jQbiS6EGK2hZfZSiI9pMdaee2Q0G7QHXpPbJjVX3kn2cAIL47Lzt2Dfxlxz6KxbJja5QqiL2ssO0uCmJPKGp5nsy2i1GFhYpavvNYimCbWIpARERE1kNJGYJ+AUsRbBVHbImIiMh66GtsiYzgiC0RERE1OzqdDocPH0ZCQkKD4lJSUqBSqZCYmCitE0IgKSkJAQEBcHZ2RnR0NA4dOmQQFx0dDZVKZbBMmDChzrbatm1bK0alUhn02dj9KpUKb775JgDg1KlTJrf57LPPGnTs1oCJLREREVmPRixF0Ol0yMjIQPfu3Q3WL1myBKmpqUhPT4dOp4NWq8WwYcNQUlJisF1cXBzOnz8vLe+//3697d28fVZWFgDgoYcekra5+f7z58/jo48+gkqlwrhx4wAAgYGBtbZ55ZVX4OrqitjYWLOP3VqwFIGIiIisx03J6Z1UWlqKyZMnY+XKlVi0aJG0XgiBtLQ0zJ8/H2PHjgUArF69Gn5+fli3bh1mzJghbevi4gKtVmt2mz4+Pga333jjDbRv3x5RUVHSulv39+WXX2LQoEFo164dAECtVtfaZtOmTRg/fjzc3NzM7ou14IgtERERNTvffPMN9u/fjylTpqC4uBjl5eV1bp+QkICRI0di6NChButPnjyJgoICxMTESOs0Gg2ioqKwb5/hrDtr166Ft7c3unbtirlz59Ya0a1LRUUFPvnkEzz++ONQqYzPWfb777/jq6++wvTp003uJycnB7m5uXVuY804YktERETWwx7K5rGtrvknMDDQYPXChQuRlJRkNGT9+vU4cOAAdDpdrfsKCgoAAH5+hlP6+fn54fTp09LtyZMnIzg4GFqtFj///DNefPFF/Pjjj1J5QX02b96MK1euYNq0aSa3Wb16Ndzd3aWRY2M+/PBDdO7cGX379jWrXWvDxJaIiIish9JSBFHzT4cOHaBWqxEXF4e4uDhoNMbnUc/Pz8fs2bORmZkJJycnk7u9dRRVCGGwLi4uTvp/aGgoQkJCEB4ejgMHDqBXr171dvvDDz9EbGwsAgICTG7z0UcfYfLkySb7ee3aNaxbtw4LFiyotz1rxcSWiIiImp2cnByzLtCQk5ODwsJChIWFSeuqqqqwa9cupKen4+jRowBqRm79/f+8qEhhYWGtUdyb9erVCw4ODjh27Fi9ie3p06fx9ddf44svvjC5ze7du3H06FFs2LDB5Daff/45ysrK8Oijj9bZnjVjjS0RERFZjzs8K8KQIUOQl5eH3NxcaQkPD8fkyZORm5uLdu3aQavVGpQUVFRUIDs7u84/9x86dAiVlZUGybApq1atgq+vL0aOHGlymw8//BBhYWHo0aNHndvcf//9tX6UZks4YktERETWQ+kFGqobtrm7uztCQ0MN1rm6usLLy0tan5iYiOTkZISEhCAkJATJyclwcXHBpEmTAAAnTpzA2rVrcd9998Hb2xuHDx/Gs88+i549e6Jfv37SfocMGYIxY8Zg5syZf3a3uhqrVq3C1KlTYW9vPG0rLi7GZ599hrffftvkcRw/fhy7du3C1q1bG/YAWBkmtkRERGQ9lNbYNjCxNce8efNw7do1xMfH4/Lly4iMjERmZibc3d0BAI6Ojvjmm2+wdOlSlJaWIjAwECNHjsTChQuhVqul/Zw4cQIXLlww2PfXX3+NM2fO4PHHHzfZ/vr16yGEwMSJE01u89FHH+Guu+4ymL3BFjGxJSIiomZHp9OZVWNrzM6dOw1uq1QqJCUlmZxVITAwENnZ2fXu99SpU7XWxcTEQAhRZ9yTTz6JJ598ss5tkpOTkZycXG8frB1rbImIiMh6NOKVx6jp44gtNSEfQ853LdUPBYpaFeEZsmO1qPsbcv3kX/VlkIJWVd1eUxANiL+Ok9/2wxuVtf3PDvLb/u24srZbb1EQvVd2ZHXpEQXtAnuRKztWfCf/9QEAN/CI7NhHMUJBy2sVxAJApILYIIVtL1EQ66qw7VEy48oB/KKwbTMprbGtslRHqCliYktERETNjpJSBGq6WIpARERE1oOlCFQHjtgSERGR9VA6KwJLEWwaE1siIiJqdliKYJtYikD12rVrF0aNGoWAgACoVCps3rzZ4H4hBJKSkhAQEABnZ2dER0fj0KFDjdNZIiKybfofj8ldHGp2w1IE28TElup19epV9OjRA+np6UbvX7JkCVJTU5Geng6dTgetVothw4ahpKTkDveUiIhsnoVqbMk2sRSB6hUbG4vY2Fij9wkhkJaWhvnz52Ps2LEAgNWrV8PPzw/r1q3DjBkzasWUl5ejvLxcul1cXHx7Ok5ERGQCSxFsE0dsSZGTJ0+ioKDA4BJ9Go0GUVFR2Ldvn9GYlJQUeHp6SktgYOCd6i4REVk7BwssYCmCrWJiS4oUFNRcHMHPz89gvZ+fn3TfrV588UUUFRVJS35+/m3vJxER2QiWIlAdWIpAFqFSqQxuCyFqrdPTaDTQaDR3oltERGRrlF55rNJSHaGmiIktKaLVagHUjNz6+/tL6wsLC2uN4hIRETUVrLG1TSxFIEWCg4Oh1WqRlZUlrauoqEB2djb69u3biD0jIiKbxCuPUR04Ykv1Ki0txfHjx6XbJ0+eRG5uLlq1aoWgoCAkJiYiOTkZISEhCAkJQXJyMlxcXDBp0qRG7DUREdkkpXWyLEWwaUxsqV4//PADBg0aJN2eM2cOAGDq1Kn4+OOPMW/ePFy7dg3x8fG4fPkyIiMjkZmZCXd398bqMhERUZ1YimCbmNhSvaKjoyGEMHm/SqVCUlISkpKS7lyniIioeVL647GKmn8iIiKgVquRkJCAhIQES/SMmgAmttRkePofB+wa/u1ZhG9T2HKa7Eg3PKmoZdWbubJjc5+T367Ie0Z+MADgquzI1f9U2LSCTzTRWtmHVynk1+K5oYeCli8riAX6YYTsWNWbyl5f4rl1smOvY6vsWCf8JDu2hpvsyGpsVtRytYKP5n5wUdT297LPlRuK2m0QpaUInO7LpjGxJSIiomaHpQi2ibMiEBERkfXgrAhUB47YEhERkfVQWmPrYKmOUFPExJaIiIiaHZYi2CaWIhAREZH1YCkC1YEjtkRERGQ9HKCsnIClCDaNiS0RERE1OyxFsE0sRSAiIiLrof/xmNzlfyO2LEWwTUxsiYiIyHpYqMZWrpSUFKhUKiQmJkrrhBBISkpCQEAAnJ2dER0djUOHDhnERUdHQ6VSGSwTJkyot71ly5YhODgYTk5OCAsLw+7duw3uT0pKQqdOneDq6oqWLVti6NCh+P777w22mTFjBtq3bw9nZ2f4+Phg9OjR+O9//yv/QWjCmNgSERGR9WjExFan0yEjIwPdu3c3WL9kyRKkpqYiPT0dOp0OWq0Ww4YNQ0lJicF2cXFxOH/+vLS8//77dba3YcMGJCYmYv78+Th48CAGDBiA2NhYnDlzRtqmY8eOSE9PR15eHvbs2YO2bdsiJiYGf/zxh7RNWFgYVq1ahSNHjmD79u0QQiAmJgZVVVXyH4wmioktERERNTs6nQ6HDx9GQoJ5l9ouLS3F5MmTsXLlSrRs2VJaL4RAWloa5s+fj7FjxyI0NBSrV69GWVkZ1q0zvKS0i4sLtFqttHh6etbZZmpqKqZPn44nnngCnTt3RlpaGgIDA7F8+XJpm0mTJmHo0KFo164dunbtitTUVBQXF+Onn/68rPSTTz6JgQMHom3btujVqxcWLVqE/Px8nDp1yqxjtyZMbImIiMh6qKuUL6gZxezUqRPefvttFBcXo7y8vM5mExISMHLkSAwdOtRg/cmTJ1FQUICYmBhpnUajQVRUFPbt22ew7dq1a+Ht7Y2uXbti7ty5tUZ0b1ZRUYGcnByD/QJATExMrf3eHJORkQFPT0/06NHD6DZXr17FqlWrEBwcjMDAwDqP2RoxsSUiIiIrUmGBBTh+/DiOHj2KuXPnwtPTEykpKSZbXL9+PQ4cOGB0m4KCAgCAn5+fwXo/Pz/pPgCYPHkyPv30U+zcuRMLFizAxo0bMXbsWJNtXrhwAVVVVfXuFwD+/e9/w83NDU5OTvj73/+OrKwseHt7G2yzbNkyuLm5wc3NDdu2bUNWVhYcHRUWHDdBnO6LiIiImp38/HyD6b40Go3J7WbPno3MzEw4OTmZ3J9KpTK4LYQwWBcXFyf9PzQ0FCEhIQgPD8eBAwfQq1cv2fsFgEGDBiE3NxcXLlzAypUr8fDDD+P777+Hr6+vtM3kyZMxbNgwnD9/Hm+99RYefvhh7N27t85jskZMbKnp6A15E2f/8wVFzf4Nh2XHvva4oqYhnstQEH1EQewlBbEAEC878hEUK2w7SHZkNUyPyJjDfYOC4FFtZIfucZEfCwCR+EJ2rHjOvPpD036UHemEJQraHa4gFgDkvzbtsLz+jeqMlz+K9p+5zyhqG28Z//N1/coBfKOsbbP9OeoqPx4YMmQI1Go1EhIS6qyzzcnJQWFhIcLCwqR1VVVV2LVrF9LT03H06FEANSO3/v7+0jaFhYW1Rltv1qtXLzg4OODYsWNGE1tvb2+o1epao7PG9uvq6ooOHTqgQ4cO6N27N0JCQvDhhx/ixRdflLbx9PSEp6cnQkJC0Lt3b7Rs2RKbNm3CxIkTTfbRGrEUgYiIiKxIJZSVIVQ2qLUhQ4YgLy8Pubm50hIeHo7JkycjNzcX7dq1g1arRVZWlhRTUVGB7Oxs9O3b1+R+Dx06hMrKSoNk+GaOjo4ICwsz2C8AZGVl1blfoGZUt76aYXO2sUYcsSUiIqJmx9wrj7m7uyM0NNRgnaurK7y8vKT1iYmJSE5ORkhICEJCQpCcnAwXFxdMmjQJAHDixAmsXbsW9913H7y9vXH48GE8++yz6NmzJ/r16yftd8iQIRgzZgxmzpwJAJgzZw6mTJmC8PBw9OnTBxkZGThz5gyeeuopADU/BHv99ddx//33w9/fHxcvXsSyZcvw22+/4aGHHgIA/Prrr9iwYQNiYmLg4+ODs2fPYvHixXB2dsZ9992n/IFsYpjYEhERkRWxTClCRESEWaUI5pg3bx6uXbuG+Ph4XL58GZGRkcjMzIS7uzuAmtHXb775BkuXLkVpaSkCAwMxcuRILFy4EGq1WtrPiRMncOHCBen2+PHjcfHiRbz66qs4f/48QkNDsXXrVrRpU1OepFar8d///herV6/GhQsX4OXlhYiICOzevRtdu3YFADg5OWH37t1IS0vD5cuX4efnh4EDB2Lfvn0GNbi2goktERERWRHLJLZK7Ny50+C2SqVCUlISkpKSjG4fGBiI7OzsevdrbF7Z+Ph4xMcb/12Dk5MTvvii7hr6gIAAbN26td62bQUTWyIiImp2zC1FIOvCH48RERGRFbHMj8ciIiLQpUsXvPfee3e4/3Q7ccSWiIiIrEglGjqzQe14slVMbImIiMiKNH6NLTVdTGyJiIio2WGNrW1ijS0RERFZEdbYkmkcsSUiIiIrUgF511+/OZ5sFRNbIiIianZYimCbWIpAREREVkRJGcKfPzxjKYJt4ogtERERWRF9ja2SeLJVKiGEaOxOUPNWXFwMT09PYEwR4NDwPwtV/bNYUft2eEBBdHtFba8RH8iOfVT1jezYvRgsOxYAdApiE5GrqG1go4LY7QrbbqUg1lF25CJsUdAucEKoZMeuUkUoahuYpyBWqyD2koJYQFniFKSoZTu3e2XHipJPFLWdq3pEVlwpitEfnigqKrptf96XPiuKPgM8XBTsqAzwfOi29pUaD0sRiIiIyIqwFIFMYykCERERWZEKKEtfOCuCLWNiS0RERM0OZ0WwTSxFICIiIivCCzSQaRyxJSIiIitSAUCtMJ5sFRNbIiIisiJMbMk0JrZERETU7LDG1jaxxpaIiIisSKUFFtbY2iqO2BIREZEVqYSyUgReecyWMbElIiKiZoelCLaJpQhERERkRXjlMTKNI7ZERERkRSqgbFyOsyLYMia2RERE1OywFME2sRSBiIiIrAivPEamccSWiIiIrEgFAJXCeLJVTGypySja9Dg84CAjsovClifIjnwITypqeaJqpezYCgVvzv0wWHZsTfwgBdFuitoG/iU7UvVbrqKWF7eWHzsPcbJjX8YI+Q0D+EIlFESfVtQ28IHsSCc8KDv2Ou6VHVvjkoJYZee4+CFXdux/VI8oarsHXpMVV4zritptDCxFsE0sRSAiIiIr0rizIqSkpEClUiExMVFaJ4RAUlISAgIC4OzsjOjoaBw6dMhovBACsbGxUKlU2Lx5c71tRUREwN3dHb6+vnjggQdw9OhRg21KS0sxc+ZMtG7dGs7OzujcuTOWL19usM2JEycwZswY+Pj4wMPDAw8//DB+//33Bh23tWBiS0RERFbEMjW2cuh0OmRkZKB79+4G65csWYLU1FSkp6dDp9NBq9Vi2LBhKCkpqbWPtLQ0qFTmlVJkZ2cjISEB+/fvR1ZWFm7cuIGYmBhcvXpV2uaZZ57Btm3b8Mknn+DIkSN45plnMGvWLHz55ZcAgKtXryImJgYqlQo7duzA3r17UVFRgVGjRqG6ulr2Y9FUsRSBiIiImp2GliKUlpZi8uTJWLlyJRYtWiStF0IgLS0N8+fPx9ixYwEAq1evhp+fH9atW4cZM2ZI2/74449ITU2FTqeDv79/vW1u27bN4PaqVavg6+uLnJwcDBw4EADw3XffYerUqYiOjgYAPPnkk3j//ffxww8/YPTo0di7dy9OnTqFgwcPSse7atUqtGrVCjt27MDQoUPNfgysAUdsiYiIyIpYphQhLCwMnTp1wttvv43i4mKUl5fX2WpCQgJGjhxZKxE8efIkCgoKEBMTI63TaDSIiorCvn37pHVlZWWYOHEi0tPTodVqZR15UVERAKBVq1bSuv79+2PLli04e/YshBD49ttv8csvv2D48OEAgPLycqhUKmg0GinGyckJdnZ22LNnj6x+NGVMbImIiMiKWCaxPX78OI4ePYq5c+fC09MTKSkpJltcv349Dhw4YHSbgoICAICfn5/Bej8/P+k+oKZkoG/fvhg9erScg4YQAnPmzEH//v0RGhoqrX/nnXfQpUsXtG7dGo6OjhgxYgSWLVuG/v37AwB69+4NV1dXPP/88ygrK8PVq1fx3HPPobq6GufPn5fVl6aMpQhERERkRW5A2XRfNwAAHTp0gFqtRlxcHOLi4gxGNG+Wn5+P2bNnIzMzE05OTib3emvdrBBCWrdlyxbs2LEDBw8elN3rmTNn4qeffqo1yvrOO+9g//792LJlC9q0aYNdu3YhPj4e/v7+GDp0KHx8fPDZZ5/hr3/9K9555x3Y2dlh4sSJ6NWrF9Rqtez+NFVMbImIiKjZycnJMavGNicnB4WFhQgLC5PWVVVVYdeuXUhPT5dmKSgoKDComy0sLJRGcXfs2IETJ06gRYsWBvseN24cBgwYgJ07d9bZh1mzZmHLli3YtWsXWrf+c97Ba9eu4aWXXsKmTZswcuRIAED37t2Rm5uLt956SyqbiImJwYkTJ3DhwgXY29ujRYsW0Gq1CA4Orvf4rQ0TWyIiIrIiFQCUzM3855XH1Go1EhISkJCQYHLrIUOGIC8vz2DdY489hk6dOuH5559Hu3btoNVqkZWVhZ49e9b0sKIC2dnZWLx4MQDghRdewBNPPGGwj27duuHvf/87Ro0aZbJtIQRmzZqFTZs2YefOnbUS0crKSlRWVsLOzrCyVK1WG53xwNvbG0BNol1YWIj777/fZNvWioktERERWRHLJLbmcnd3N6hpBQBXV1d4eXlJ6xMTE5GcnIyQkBCEhIQgOTkZLi4umDRpEgBAq9Ua/cFYUFCQQbI6ZMgQjBkzBjNnzgRQ84O1devW4csvv4S7u7tUs+vp6QlnZ2d4eHggKioKzz33HJydndGmTRtkZ2djzZo1SE1Nlfa7atUqdO7cGT4+Pvjuu+8we/ZsPPPMM7j77rsb9FhYAya2RERE1OxY8spj8+bNw7Vr1xAfH4/Lly8jMjISmZmZcHd3b9B+9OUCevoLLein8tJbtWoVpk2bBqDmh20vvvgiJk+ejEuXLqFNmzZ4/fXX8dRTT0nbHz16FC+++CIuXbqEtm3bYv78+XjmmWfkHWwTx8SWiIiIrEgllI3Y1vx4zNxSBGNurYlVqVRISkpCUlKS2fsQovYxnDp1qt5tbqXVarFq1ao6t3njjTfwxhtvmN03a8bEloiIiKxIBQAlV8y6YamOUBPExJaIiIiaHUuWIlDTwQs0EBERkRWxzAUaIiIi0KVLF7z33nt3uP90O3HElpqQZAANK7Sv8a3Cdn+XHfkZFihse6+C2EsKYjMUxAKAkililLWtEtmyY48rmdMdQHuMUxDdQXbkMqxU0C4Qj8UKotsrahvQyY4s/l8CIk+iglgAWCI7UiW6KWr5fgXnaYSicxSQ/56kpDSgoSoVtldlqY5QE8TEloiIiJodliLYJpYiEBERkRVhKQKZxhFbIiIisiIVANQK4lmKYMuY2BIREZEVqYSy5PRO1gPTncZSBKrXrl27MGrUKAQEBEClUmHz5s0G90+bNg0qlcpg6d27d+N0loiIyAw6nQ6HDx9u8MUZqGljYkv1unr1Knr06IH09HST24wYMQLnz5+Xlq1bt97BHhIRUfPBGlsyjaUIVK/Y2FjExsbWuY1Go4FWqzVrf+Xl5SgvL5duFxcXK+ofERE1J5VQNi7HUgRbxsSWLGLnzp3w9fVFixYtEBUVhddffx2+vr5Gt01JScErr7xyh3tIRET0J073ZZtYikCKxcbGYu3atdixYwfefvtt6HQ6DB482GBU9mYvvvgiioqKpCU/P/8O95iIiKxXJZSVIVQCYCmCreKILSk2fvx46f+hoaEIDw9HmzZt8NVXX2Hs2LG1ttdoNNBoNHeyi0REZCNUuAEV5F+eTUBAWLA/1LQwsSWL8/f3R5s2bXDs2LHG7goREZFRLEWwTSxFIIu7ePEi8vPz4e/v39hdISIiG2OHG4oXgKUItoojtlSv0tJSHD9+XLp98uRJ5ObmolWrVmjVqhWSkpIwbtw4+Pv749SpU3jppZfg7e2NMWPGNGKviYjIFtlZoBSB1x6zXUxsqV4//PADBg0aJN2eM2cOAGDq1KlYvnw58vLysGbNGly5cgX+/v4YNGgQNmzYAHd398bqMhERUZ1YimCbmNhSvaKjoyGE6VL77du338HeEBFRc2aHKgXjtZB+OBYREQG1Wo2EhARefcyGMLGlJmQcALWMuPYK2z2iIHZQ/ZvUQbV0h+xYMXub7NhQ3Cs7FgB+hnkX4zBumaK2hapUQfS3itq+js9kxzrhA9mxj8qO1PuP7MgjeF5Ry8GoPTOKuZywU0HLbgpigRnIlR27TUnWBWA4+iuILlDWuOzXdhWAPxS2bR5HQHFia3wySrIFTGyJiIio2WEpgm3irAhERERkNRxQM2ord3H43344K4Jt4ogtERERWQ1HKBuVq7ZUR6hJYmJLREREVoOJLdWFiS0RERE1O6yxtU2ssSUiIiKrYY+aOlm5i35EjzW2tokjtkRERGQ1HCFvYkg9XnXMtnHEloiIiJodnU6Hw4cPN/jiDCkpKVCpVEhMTJTWCSGQlJSEgIAAODs7Izo6GocOHTKImzFjBtq3bw9nZ2f4+Phg9OjR+O9//1tnW0lJSVCpVAaLVmt6ruEZM2ZApVIhLS3NYH10dHSt/UyYMKFBx20tmNgSERGR1VAy1Zd+AeSVIuh0OmRkZKB79+4G65csWYLU1FSkp6dDp9NBq9Vi2LBhKCkpkbYJCwvDqlWrcOTIEWzfvh1CCMTExKCqqu4x5K5du+L8+fPSkpeXZ3S7zZs34/vvv0dAQIDR++Pi4gz28/7775t93NaEiS0RERFZDUvNY9tQpaWlmDx5MlauXImWLVtK64UQSEtLw/z58zF27FiEhoZi9erVKCsrw7p166TtnnzySQwcOBBt27ZFr169sGjRIuTn5+PUqVN1tmtvbw+tVistPj4+tbY5e/YsZs6cibVr18LBwfgRuri4GOzH09NT3gPRxDGxJSIiombnm2++wf79+zFlyhQUFxejvLzuC+0mJCRg5MiRGDp0qMH6kydPoqCgADExMdI6jUaDqKgo7Nu3z+i+rl69ilWrViE4OBiBgYF1tnvs2DEEBAQgODgYEyZMwK+//mpwf3V1NaZMmYLnnnsOXbt2NbmftWvXwtvbG127dsXcuXMNRpNtCX88RkRERFbDEcqSF/2I3q0J5cKFC5GUlGQ0Zv369Thw4AB0Ol2t+woKCgAAfn5+Buv9/Pxw+vRpg3XLli3DvHnzcPXqVXTq1AlZWVlwdHSEKZGRkVizZg06duyI33//HYsWLULfvn1x6NAheHl5AQAWL14Me3t7PP300yb3M3nyZAQHB0Or1eLnn3/Giy++iB9//BFZWVkmY6wVE1siIiKyGpZKbDt06AC1Wo24uDjExcVBo9EY3T4/Px+zZ89GZmYmnJycTO5XpVIZ3BZC1Fo3efJkDBs2DOfPn8dbb72Fhx9+GHv37jW539jYWOn/3bp1Q58+fdC+fXusXr0ac+bMQU5ODpYuXYoDBw7UautmcXFx0v9DQ0MREhKC8PBwHDhwAL169TIZZ42Y2BIREVGzk5OTY9YFGnJyclBYWIiwsDBpXVVVFXbt2oX09HQcPXoUQM3Irb+/v7RNYWFhrVFcT09PeHp6IiQkBL1790bLli2xadMmTJw40aw+u7q6olu3bjh27BgAYPfu3SgsLERQUJBB35599lmkpaWZrN/t1asXHBwccOzYMSa2RLdN0V5AzlVgVFsUNvwvBbHDFbV8bbaSiztWyI50VdBqjVGyI1VYpKhlL2yTHdsKkxS1/Qs+VhB9RHakG+5R0C7gjVzZsb8pvACpE4bWv5FJSqYj+l1BLPA+liiITlfU9hfYIzv2ARQratsOdU8/ZVopgCGK2jaX/kILcunHNSMiIqBWq5GQkFDnlF9DhgypNRPBY489hk6dOuH5559Hu3btoNVqkZWVhZ49ewIAKioqkJ2djcWLF9fZFyFEvbW9NysvL8eRI0cwYMAAAMCUKVNq1fwOHz4cU6ZMwWOPPWZyP4cOHUJlZaVBIm4rmNgSERGR1VAyswHwZ2JrLnd3d4SGhhqsc3V1hZeXl7Q+MTERycnJCAkJQUhICJKTk+Hi4oJJk2q+TP/666/YsGEDYmJi4OPjg7Nnz2Lx4sVwdnbGfffdJ+13yJAhGDNmDGbOnAkAmDt3LkaNGoWgoCAUFhZi0aJFKC4uxtSpUwEAXl5eUq2tnoODA7RaLe6++24AwIkTJ7B27Vrcd9998Pb2xuHDh/Hss8+iZ8+e6NevXwMfjaaPiS0RERFZjTud2Jpj3rx5uHbtGuLj43H58mVERkYiMzMT7u7uAAAnJyfs3r0baWlpuHz5Mvz8/DBw4EDs27cPvr6+0n5OnDiBCxcuSLd/++03TJw4ERcuXICPjw969+6N/fv3o02bNmb3zdHREd988w2WLl2K0tJSBAYGYuTIkVi4cCHUaiXXcGuamNgSERFRs6PT6cyqsTVm586dBrdVKhWSkpJMzqoQEBCArVu31rvfW2ti169f3+C+3bqPwMBAZGdnN3g/1orz2BIREZHVaMwrj1HTxxFbIiIishr6K48RGcPEloiIiJodJaUI1HSxFIGIiIishoMFFoClCLaKI7ZERERkNW6ukyW6FRNbIiIianZYimCbWIpAREREVkP/4zG5C0sRbBtHbImIiMhqKC1FEJbqCDVJTGyJiIio2WEpgm1iKQIRERFZDV6ggerCEVsiIiKyGkov0FBtqY5Qk8TElpoOzwgAahmBgxQ2rJUduR73K2p5Ap5XEN1KduSPCvsN7JUdKZCosO11CmLl9xsAVD7yR3bEHzsUtHxJQSxwTEGsEx5Q1PZ1fK2g7X2yY1Wim+xYAHhStUB2bEayoqYhXhosO1Yl7lPWtuqozMgKRe02BpYi2CaWIhAREZHVYCkC1YUjtkRERGQ1lM6KwFIE28bEloiIiKyG0hrbKkt1hJokJrZERETU7LDG1jaxxpaIiIisBmtsqS4csSUiIiKrobTGlqUIto2JLRERETU7LEWwTSxFICIiIqthj5ofkMld9CN6LEWwTRyxJSIiIquhtBThhqU6Qk0SE1siIiJqdliKYJtYikBERERWg7MiUF04YktERERWQ+kFGiot1RFqkpjYEhERUbPDUgTbxFIEIiIishosRaC6cMSWiIiIrIbSWRFYimDbmNhSk3EUOrih4X8W6oxnFLY8XHbkgyhT1LLKc7HsWFE0UXbsddwrO7aGVkHsPEUtz4G/7NiJmKSobfHHDtmx32Kw7NhB6Cw7FoCCRwy4jsOK2rbHdQXR22VHCtU4Be0CSs7x919KUNTyAMg/z7BVUdOoRrGCuDXKGjeT0hrbCoXtp6Sk4KWXXsLs2bORlpYGABBC4JVXXkFGRgYuX76MyMhIvPfee+jatSsA4NKlS1i4cCEyMzORn58Pb29vPPDAA3jttdfg6elpsq0bN24gKSkJa9euRUFBAfz9/TFt2jS8/PLLsLOzM6ttAMjIyMC6detw4MABlJSU4PLly2jRooXCR6JpYikCERERNTs6nQ6HDx9GQoL5X0R0Oh0yMjLQvXt3g/VLlixBamoq0tPTodPpoNVqMWzYMJSUlAAAzp07h3PnzuGtt95CXl4ePv74Y2zbtg3Tp0+vs73FixdjxYoVSE9Px5EjR7BkyRK8+eabePfdd81uGwDKysowYsQIvPTSS2Yfq7XiiC0RERFZDaWlCPrYsLAwqNVqxMXFIS4uDhqNBhqNxmRcaWkpJk+ejJUrV2LRokXSeiEE0tLSMH/+fIwdOxYAsHr1avj5+WHdunWYMWMGQkNDsXHjRimmffv2eP311/HII4/gxo0bsLc3no599913GD16NEaOHAkAaNu2LT799FP88MMPZrcNAImJiQCAnTt3NvjxsjYcsSUiIiKrYakfjx0/fhxHjx7F3Llz4enpiZSUlDrbTUhIwMiRIzF06FCD9SdPnkRBQQFiYmKkdRqNBlFRUdi3b5/J/RUVFcHDw8NkUgsA/fv3xzfffINffvkFAPDjjz9iz549uO+++xS1bcs4YktERETNTn5+vsF0X3WN1q5fvx4HDhyATqerdV9BQQEAwM/Pz2C9n58fTp8+bXR/Fy9exGuvvSaNqJry/PPPo6ioCJ06dYJarUZVVRVef/11TJw4UXbbto6JLREREVkNpT8ec/jfv0OGDIFarUZCQkKddbb5+fmYPXs2MjMz4eTkZHI7lUplcFsIUWsdABQXF2PkyJHo0qULFi5cWGdfN2zYgE8++QTr1q1D165dkZubi8TERAQEBGDq1KkNbrs5YGJLREREVsMBfyancuMbIicnB4WFhQgLC5PWVVVVYdeuXUhPT8fRo0cBQJq1QK+wsLDWSGpJSQlGjBgBNzc3bNq0CQ4OdffmueeewwsvvIAJEyYAALp164bTp08jJSUFU6dOhVarNbvt5oI1tkRERNTsmDsrwpAhQ5CXl4fc3FxpCQ8Px+TJk5Gbm4t27dpBq9UiKytLiqmoqEB2djb69u0rrSsuLkZMTAwcHR2xZcuWOkd/9crKyqRpvfTUajWqq6sBAMHBwWa13ZxwxJaIiIishqVmRYiIiDCrFMHd3R2hoaEG61xdXeHl5SWtT0xMRHJyMkJCQhASEoLk5GS4uLhg0qSaubNLSkoQExODsrIyfPLJJyguLkZxcc2cwT4+PlCr1QBqkugxY8Zg5syZAIBRo0bh9ddfR1BQELp27YqDBw8iNTUVjz/+OICaEoT62gZqRnQLCgpw/PhxAEBeXh7c3d0RFBSEVq1aKXg0mx4mtkRERGQ1LFVja0nz5s3DtWvXEB8fL10kITMzE+7u7gBqyhm+//57AECHDh0MYk+ePIm2bdsCAE6cOIELFy5I97377rtYsGAB4uPjUVhYiICAAMyYMQN/+9vfzG4bAFasWIFXXnlFuj1w4EAAwKpVqzBt2jSLPhaNTSWEEI3dCWreiouL4enpif+gqJGuPCb/CkU30EtRyw6eLrJjlVx5DDihIBZQduWx9xW1rOzKY8pEKLgilLIrj52XHQsATgoes+voUP9GdbiBn2TH2qPu6ZfqpvTKYyMUxCq98tgC2bF7vlLUNKpGyrvyWDGK0RKB0hRWt4P+s+JUEaCkieJioK0nbmtfqfGwxpaIiIishqXmsY2IiECXLl3w3nvv3dH+0+3FUgQiIiKyGvbVNYuSeLJdTGyJiIio2dHpdCxFsEEsRSAiIiKrYXdD+QKwFMFWccSWmoy7MQseMn7rqgpbqahdkXOP7Fh7dFbUdkHRpwqih8uOVM2dpqBdQLwVJzt2jYIfMgFAKv6jIPpJRW1/gVzZsWOxR0HLyiZa/0xRtJLfnwP2Cn4s+Bhekx27SsGP9WqMUhA7SFHLu5X8IHaksudL7vuKHa4qbLcBbd2UnMqNJ9vFxJaIiIisBhNbqgsTWyIiImp2WGNrm1hjS0RERFZDdQNQVSpYWGNr0zhiS0RERNaj4n+LkniyWRyxpTqlpKQgIiIC7u7u8PX1xQMPPICjR48abCOEQFJSEgICAuDs7Izo6GgcOnSokXpMRERUP51Oh8OHDyMhQdmV4qhpYWJLdcrOzkZCQgL279+PrKws3LhxAzExMbh69c9fwC5ZsgSpqalIT0+HTqeDVqvFsGHDUFJS0og9JyIim1RhgQUsRbBVLEWgOm3bts3g9qpVq+Dr64ucnBwMHDgQQgikpaVh/vz5GDt2LABg9erV8PPzw7p16zBjxoxa+ywvL0d5ebl0u7hY3rXJiYioGar836IknmwWR2ypQYqKigAArVq1AgCcPHkSBQUFiImJkbbRaDSIiorCvn37jO4jJSUFnp6e0hIYGHj7O05ERHQTliLYJia2ZDYhBObMmYP+/fsjNDQUAFBQUAAA8PMznEDez89Puu9WL774IoqKiqQlPz//9naciIhsRyWUlSH8b8SWpQi2iaUIZLaZM2fip59+wp49ta+gpFKpDG4LIWqt09NoNNBoNLelj0REZOM4KwLVgYktmWXWrFnYsmULdu3ahdatW0vrtVotgJqRW3//Py+VWlhYWGsUl4iIqKngBRpsE0sRqE5CCMycORNffPEFduzYgeDgYIP7g4ODodVqkZWVJa2rqKhAdnY2+vbte6e7S0REto6zIlAdOGJLdUpISMC6devw5Zdfwt3dXaqb9fT0hLOzM1QqFRITE5GcnIyQkBCEhIQgOTkZLi4umDRpUiP3noiIbI6+xlZJPNksJrZUp+XLlwMAoqOjDdavWrUK06ZNAwDMmzcP165dQ3x8PC5fvozIyEhkZmbC3d39DveWiIhsXgUAB4XxZLOY2FKdhBD1bqNSqZCUlISkpKTb3yEiIiILYI2tbWJiS01GW7wLOzT8TUbkbKt/ozqNkx2pKlygqGXhO0J+20/JP27x1mLZsQCgenil7Nh7/6moaTyKaQqipypqu72CWNWv/WXHTmqnoGEAa6HkNRKvrHGslx25CvJ/gPoBdsiOBYAnoOS1/aCitgE3BbHKznHgW5lx5fVvYikWGrGNiIiAWq1GQkIC57K1IUxsiYiIyHqwxpbqwMSWiIiImh2WItgmTvdFRERE1oPTfVEdOGJLRERE1qMCyrIXzopg05jYEhERUbPDUgTbxFIEIiIish76H4/JXf734zG5pQgpKSnSxYn0hBBISkpCQEAAnJ2dER0djUOHDkn3X7p0CbNmzcLdd98NFxcXBAUF4emnn0ZRUVG97S1btgzBwcFwcnJCWFgYdu/ebXD/tGnToFKpDJbevXtL9586darW/frls88+a9CxWwMmtkRERGQ9Ki2wyKTT6ZCRkYHu3bsbrF+yZAlSU1ORnp4OnU4HrVaLYcOGoaSkBABw7tw5nDt3Dm+99Rby8vLw8ccfY9u2bZg+fXqd7W3YsAGJiYmYP38+Dh48iAEDBiA2NhZnzpwx2G7EiBE4f/68tGzdulW6LzAw0OC+8+fP45VXXoGrqytiY2PlPxhNFEsRiIiIqNlpaClCaWkpJk+ejJUrV2LRokXSeiEE0tLSMH/+fIwdOxYAsHr1avj5+WHdunWYMWMGQkNDsXHjRimmffv2eP311/HII4/gxo0bsLc3no6lpqZi+vTpeOKJJwAAaWlp2L59O5YvX46UlBRpO41GA61Wa3QfarW61n2bNm3C+PHj4eamZM7kpokjtkRERGQ9LDQrQlhYGDp16oS3334bxcXFKC+v+yITCQkJGDlyJIYOHWqw/uTJkygoKEBMTIy0TqPRICoqCvv27TO5v6KiInh4eJhMaisqKpCTk2OwXwCIiYmptd+dO3fC19cXHTt2RFxcHAoLC022m5OTg9zc3HpHi60VE1siIiKyHhaqsT1+/DiOHj2KuXPnwtPT02AE9Fbr16/HgQMHjG5TUFAAAPDzM7xSnp+fn3TfrS5evIjXXnsNM2bMMNnmhQsXUFVVVe9+Y2NjsXbtWuzYsQNvv/02dDodBg8ebDJR//DDD9G5c2f07dvXZNvWjKUIREREZD0qoGxY7n8jth06dIBarUZcXBzi4uKg0WiMbp6fn4/Zs2cjMzMTTk5OJnerUqkMbgshaq0DgOLiYowcORJdunTBwoUL6+1uffsdP3689P/Q0FCEh4ejTZs2+Oqrr6TSCL1r165h3bp1WLBA2eXgmzImtkRERNTs5OTkmFVjm5OTg8LCQoSFhUnrqqqqsGvXLqSnp+Po0aMAakZu/f39pW0KCwtrjbaWlJRgxIgRcHNzw6ZNm+Dg4GCyXW9vb6jV6lqjvsb2ezN/f3+0adMGx44dq3Xf559/jrKyMjz66KN1H7QVYykCERERWY87fOWxIUOGIC8vD7m5udISHh6OyZMnIzc3F+3atYNWq0VWVtafXayoQHZ2tsGf+4uLixETEwNHR0ds2bKlztFfAHB0dERYWJjBfgEgKyurzjKCixcvIj8/3yDJ1vvwww9x//33w8fHp862rRlHbImIiMh6VELZsFwDp/tyd3dHaGiowTpXV1d4eXlJ6xMTE5GcnIyQkBCEhIQgOTkZLi4umDRpEoCakdqYmBiUlZXhk08+QXFxMYqLiwEAPj4+UKvVAGqS6DFjxmDmzJkAgDlz5mDKlCkIDw9Hnz59kJGRgTNnzuCpp54CUDNTQ1JSEsaNGwd/f3+cOnUKL730Ery9vTFmzBiDPh8/fhy7du0ymArMFjGxpSbjFB6HB0z/WcYUVbdPFbUr8trLD76kqGkAE2RHihWDZcd+gB2yYwFg4z/lx8bUv0k9XpYduQaTFLX8KGrXy5lLtJN/nl3HL7JjAWCJkD9X5TxVvKK2gXkKYnWyIzfWv0mdnoCS4z6sqG2VOCg7Vqj2Kmz7RXmBxVcBzzcVtX2nWfLKY/PmzcO1a9cQHx+Py5cvIzIyEpmZmXB3dwdQU87w/fffA6ip7b3ZyZMn0bZtWwDAiRMncOHCBem+8ePH4+LFi3j11Vdx/vx5hIaGYuvWrWjTpg2Amqm88vLysGbNGly5cgX+/v4YNGgQNmzYILWt99FHH+Guu+6qNcuCrWFiS0RERNajAlDwHdOgFEGtViMhIQEJCQkN2sXOnTsNbqtUKiQlJSEpKcno9tHR0RBC1LvfU6dO1VoXHx+P+HjjX7ScnZ2xffv2evcLAMnJyUhOTjZrW2vGxJaIiIish4USW7JNTGyJiIio2bFkKQI1HZwVgYiIiKyHhS7QYO6sCGRdOGJLRERE1kNpKQFLEWwaE1siIiJqdliKYJtYikBERETWo9ICC1iKYKs4YktERETWo4EXWLB4PDVpTGyJiIio2WEpgm1iKQIRERFZDyUzIugXsBTBVnHEloiIiKxHBYD6L+JlGksRbBoTWyIiIrIelVCW2N6wVEeoKWJiS0RERM0Oa2xtE2tsiYiIyHqwxpbqwBFbIiIish4VAKoVxLMUwaYxsaUmoxAf4Roa/mchkXdRUbt/Q4jsWNFpqqK2AZ3syA9wWHbsE+giO7bGOAWx2xW2/arsyEfRX2HbfrIjVeJp2bFClSk7FgD2qpQUJD6jqO3r+IvsWCc8KTt228jXZMcCAL56THbot9imqOkqVbGC6BcVtS1UebLiinENnopavvNYimCbWIpARERE1qMSysoQeOUxm8YRWyIiIrIeFQDUCuKrLNURaoqY2BIREVGzw1IE28RSBCIiIrIenBWB6sARWyIiIrIelVBWTqBkRgVq8pjYEhERUbPDUgTbxFIEIiIish6VFljAUgRbxRFbIiIish4VUDYsx1IEm8bEloiIiKxHJQCVgngl1yuhJo+JLRERETU7rLG1TayxJSIiIuvRyNN9paSkQKVSITExUVonhEBSUhICAgLg7OyM6OhoHDp0yCAuIyMD0dHR8PDwgEqlwpUrV+pta/ny5ejevTs8PDzg4eGBPn364P/+7/9Mbj9jxgyoVCqkpaVJ6y5duoRZs2bh7rvvhouLC4KCgvD000+jqKioQcdtLZjYEhERkfWwUGIrh06nQ0ZGBrp3726wfsmSJUhNTUV6ejp0Oh20Wi2GDRuGkpISaZuysjKMGDECL730ktnttW7dGm+88QZ++OEH/PDDDxg8eDBGjx5dK2kGgM2bN+P7779HQECAwfpz587h3LlzeOutt5CXl4ePP/4Y27Ztw/Tp0xt49NaBpQhERETU7DS0FKG0tBSTJ0/GypUrsWjRImm9EAJpaWmYP38+xo4dCwBYvXo1/Pz8sG7dOsyYMQMApBHenTt3mt3mqFGjDG6//vrrWL58Ofbv34+uXbtK68+ePYuZM2di+/btGDlypEFMaGgoNm7cKN1u3749Xn/9dTzyyCO4ceMG7O1tKxW0raMhqyRETSV/CYplxTujpP6N6lAOB9mxxUq++gNQMsv4NZmPFwAUK75YermCWKVtX1UQe0Nh2wp+Tl18XX6oomMGKhWdK0qea+C6grYrlDzelfLbBYBiBefKVQXHXNO2/Hg7xe9J12RFFaPm/Na/n99OxZUK4//3b1hYGNRqNeLi4hAXFweNRgONRmMyLiEhASNHjsTQoUMNEtuTJ0+ioKAAMTEx0jqNRoOoqCjs27dPSmyVqqqqwmeffYarV6+iT58+0vrq6mpMmTIFzz33nEGyW5eioiJ4eHjYXFILMLGlJkD/p5peCGzknjTckkZt3VN25NOK235T8R7ke7AR21bA83n5oRbsxp1ve7kFeiFDprKeN2a0NSspKYGn5+05fkdHR2i1WgQWFCjel5ubG44fPw4AmDt3LubOnYuFCxciKSnJ6Pbr16/HgQMHoNPpat1X8L/++Pn5Gaz38/PD6dOnFfc1Ly8Pffr0wfXr1+Hm5oZNmzahS5cu0v2LFy+Gvb09nn7avHf2ixcv4rXXXrNYwt3UMLGlRhcQEID8/Hy4u7tDpao9h0txcTECAwORn5/PX7CaiY9Zw/Exazg+Zg1nq4+ZEAIlJSW16jstycnJCSdPnkRFhdJR6Zr+3vp5Y2q0Nj8/H7Nnz0ZmZiacnJxM7vPW/RlrQ467774bubm5uHLlCjZu3IipU6ciOzsbXbp0QU5ODpYuXYoDBw6Y1VZxcTFGjhyJLl26YOHChYr71hQxsaVGZ2dnh9atW9e7nf5XoWQ+PmYNx8es4fiYNZwtPma3a6T2Zk5OTnUml7dDTk4OCgsLERYWJq2rqqrCrl27kJ6ejqNHjwKoGbn19/eXtiksLKw1iiuHo6MjOnToAAAIDw+HTqfD0qVL8f7772P37t0oLCxEUFCQQd+effZZpKWl4dSpU9L6kpISjBgxQhr1dXCQX4bXlDGxJSIiIjJhyJAhyMvLM1j32GOPoVOnTnj++efRrl07aLVaZGVloWfPngCAiooKZGdnY/HixRbvjxAC5eU1te9TpkzB0KFDDe4fPnw4pkyZgscee0xaV1xcjOHDh0Oj0WDLli13/MvBncTEloiIiMgEd3d3hIaGGqxzdXWFl5eXtD4xMRHJyckICQlBSEgIkpOT4eLigkmTJkkxBQUFKCgokGp78/Ly4O7ujqCgILRq1QpATRI9ZswYzJw5EwDw0ksvITY2FoGBgSgpKcH69euxc+dObNu2DQDg5eUFLy8vg745ODhAq9Xi7rvvBlAzUhsTE4OysjJ88sknKC4uRnFxzU/ofHx8oFarLf2QNSomttTkaTQaLFy4sM5fq5IhPmYNx8es4fiYNRwfM9s0b948XLt2DfHx8bh8+TIiIyORmZkJd3d3aZsVK1bglVdekW4PHDgQALBq1SpMmzYNAHDixAlcuHBB2ub333/HlClTcP78eXh6eqJ79+7Ytm0bhg0bZnbfcnJy8P333wOAVNKgd/LkSbRt27ahh9ukqcSdmJuDiIiIiOg245XHiIiIiMgmMLElIiIiIpvAxJaIiIiIbAITWyIiIiKyCUxsiYiIiMgmMLGlJm/ZsmUIDg6Gk5MTwsLCsHv37sbuUpOVlJQElUplsGi12sbuVpOya9cujBo1CgEBAVCpVNi8ebPB/UIIJCUlISAgAM7OzoiOjsahQ4cap7NNRH2P2bRp02qdd717926czjYBKSkpiIiIgLu7O3x9ffHAAw9IV6fS43lGdHswsaUmbcOGDUhMTMT8+fNx8OBBDBgwALGxsThz5kxjd63J6tq1K86fPy8tt14xp7m7evUqevTogfT0dKP3L1myBKmpqUhPT4dOp4NWq8WwYcNQUlJyh3vadNT3mAHAiBEjDM67rVu33sEeNi3Z2dlISEjA/v37kZWVhRs3biAmJgZXr16VtuF5RnSbCKIm7N577xVPPfWUwbpOnTqJF154oZF61LQtXLhQ9OjRo7G7YTUAiE2bNkm3q6urhVarFW+88Ya07vr168LT01OsWLGiEXrY9Nz6mAkhxNSpU8Xo0aMbpT/WoLCwUAAQ2dnZQgieZ0S3E0dsqcmqqKhATk4OYmJiDNbHxMRg3759jdSrpu/YsWMICAhAcHAwJkyYgF9//bWxu2Q1Tp48iYKCAoNzTqPRICoqiudcPXbu3AlfX1907NgRcXFxKCwsbOwuNRlFRUUAIF02lecZ0e3DxJaarAsXLqCqqgp+fn4G6/38/FBQUNBIvWraIiMjsWbNGmzfvh0rV65EQUEB+vbti4sXLzZ216yC/rziOdcwsbGxWLt2LXbs2IG3334bOp0OgwcPRnl5eWN3rdEJITBnzhz0798foaGhAHieEd1O9o3dAaL6qFQqg9tCiFrrqEZsbKz0/27duqFPnz5o3749Vq9ejTlz5jRiz6wLz7mGGT9+vPT/0NBQhIeHo02bNvjqq68wduzYRuxZ45s5cyZ++ukn7Nmzp9Z9PM+ILI8jttRkeXt7Q61W1xrBKCwsrDXSQca5urqiW7duOHbsWGN3xSroZ5DgOaeMv78/2rRp0+zPu1mzZmHLli349ttv0bp1a2k9zzOi24eJLTVZjo6OCAsLQ1ZWlsH6rKws9O3bt5F6ZV3Ky8tx5MgR+Pv7N3ZXrEJwcDC0Wq3BOVdRUYHs7Gyecw1w8eJF5OfnN9vzTgiBmTNn4osvvsCOHTsQHBxscD/PM6Lbh6UI1KTNmTMHU6ZMQXh4OPr06YOMjAycOXMGTz31VGN3rUmaO3cuRo0ahaCgIBQWFmLRokUoLi7G1KlTG7trTUZpaSmOHz8u3T558iRyc3PRqlUrBAUFITExEcnJyQgJCUFISAiSk5Ph4uKCSZMmNWKvG1ddj1mrVq2QlJSEcePGwd/fH6dOncJLL70Eb29vjBkzphF73XgSEhKwbt06fPnll3B3d5dGZj09PeHs7AyVSsXzjOh2adQ5GYjM8N5774k2bdoIR0dH0atXL2nKHKpt/Pjxwt/fXzg4OIiAgAAxduxYcejQocbuVpPy7bffCgC1lqlTpwohaqZiWrhwodBqtUKj0YiBAweKvLy8xu10I6vrMSsrKxMxMTHCx8dHODg4iKCgIDF16lRx5syZxu52ozH2WAEQq1atkrbheUZ0e6iEEOLOp9NERERERJbFGlsiIiIisglMbImIiIjIJjCxJSIiIiKbwMSWiIiIiGwCE1siIiIisglMbImIiIjIJjCxJSIiIiKbwMSWiIiIiGwCE1siIiIisglMbImIiIjIJjCxJSIiIiKb8P9P/ll6MHg9HAAAAABJRU5ErkJggg==",
      "text/plain": [
       "<Figure size 640x480 with 2 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "plot_log1p_dlog(L,dlog_P,path='plots/eigenvectors/',title=f\"log1p of dlog of eigenvectors of uDFT for n={n} q={q} deg={L.degree()}\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "370a773c-5d21-49d0-950b-887fc5c91332",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "SageMath 10.5",
   "language": "sage",
   "name": "sagemath-10.5"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
