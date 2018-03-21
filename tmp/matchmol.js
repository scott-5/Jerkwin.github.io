var $=function(id){return document.getElementById(id)}
function trim(s){return s.replace(/^(\s|\u00A0)+/,'').replace(/(\s|\u00A0)+$/,'')}
function norm(s){return s.replace(/^\s*\n*/,'').replace(/\s*\n*$/,'').replace(/\s+[\n|$]/g,'\n')}
function printf() {
	var map = {
		s: function(str, fmt) { var n=str.length; return n>fmt?str:Array(fmt-n+1).join(' ')+str},
		f: function(str, fmt) { fmt=fmt.split('.'); str=parseFloat(str).toFixed(fmt[1]);
			var m=fmt[0], n=str.length; return n>m?str:Array(m-n+1).join(' ')+str
		}
	}
	var args = Array.prototype.slice.call(arguments).slice();
	return args.shift().toString().replace(/%(-*\d*\.*\d*)([sf])/g, function(_, fmt, type){
		if(!args.length) return 'ERROR'
		return map[type](args.shift().toString(), fmt);
	});
}

var NMAX=10, drawMode=1, G=[], P=[], adjG=[], adjP=[],
	idxFit=[], fitAdj=[], fitX=[], fitY=[], fitZ=[], isMatch=[]

var img = []

img[1]=new Image();  img[1].src = "H.png";
img[6]=new Image();  img[6].src = "C.png";
img[7]=new Image();  img[7].src = "N.png";
img[8]=new Image();  img[8].src = "O.png";
imgBond=new Image(); imgBond.src='bond.png'

//   function initialize()
//
//   Initialization Routines
//     - Load properties of elements
//     - Define Handlers for mouse events
//     - Draw molecule
//
function initialize() {

  // Define handler for drawing area (canvas)
  // activeWin defines routines to handle mouse events
  var activeCanvas = activeWin();
  var molfile = $('fileRef');

  // If no canvas element found, see if userDefined() creates one.
  if ( ! activeCanvas ) {
    userDefined();
    activeCanvas = activeWin("");
    if ( ! activeCanvas ) {
      alert("No canvas present. Either create a canvas or remove onLoad=initialize() statement from <body>.");
      }
    return;
    }

  // Handler for loading files. After file selected, calls "MyFileReader"
  if ( molfile ) {
    // Verify browser supports HTML5 FileReader
    if (window.File && window.FileReader && window.FileList && window.Blob) {
        molfile.addEventListener('change', MyFileReader, false);
      molfile = $('fileFit');
        molfile.addEventListener('change', MyFileReader, false);
      } else {
        alert('The File APIs are not fully supported by your browser.');
      }
    }

  // Create periodic table
  //drawPeriodic();

  // Initialize molecular geometry
  var molecule = Mol();

  // Set display mode to View
  viewmode();

  // Draw molecule
  if ( typeof molecule == 'undefined' )
    delMolecule();
  if (molecule[0].numatoms < 1) {
    delMolecule();
    loadMolecule();
    formula();  // Display molecular formula
    centerMolecule();  // Center molecule
    }
   drawMolecule();

  // Set default drawing atom
  //pickClouds(4);
  //pickElem("C");

  // Execute user-defined routine
  userDefined();

  }

  function userDefined() {
    activeWin("molFit",2);
    delMolecule();

  // Add atomic coordinates to molecule array
  addAtom(1,  0.874,  0.618,  0.000);
  addAtom(6,  0.000,  0.000,  0.000);
  addAtom(1, -0.874,  0.618,  0.000);
  addAtom(1,  0.000, -0.618,  0.874);
  addAtom(1,  0.000, -0.618, -0.874);
  // Add bonds to bond array
  addBond(2, 1);
  addBond(2, 3);
  addBond(2, 4);
  addBond(2, 5);
    formula();  // Display molecular formula
    centerMolecule();
    drawMolecule();
}

function match() {
	var i, j, k, mol, n, adj

	activeWin('molRef');
	mol=Mol(); n=mol[0].numatoms; adj=BondMatrix()
	for(i=0; i<n; i++) {
		G[i]=[]; adjG[i]=[]

		adjG[i][0]=mol[i+1].atomicnumber
		for(j=1; j<NMAX; j++) adjG[i][j]=0

		for(j=0; j<n; j++) {
			k=adj[i+1][j+1]
			G[i][j]=k
			adjG[i][mol[j+1].atomicnumber] += k
		}
	}

	activeWin("molFit",2);
	mol=Mol(); n=mol[0].numatoms; adj=BondMatrix()
	for(i=0; i<n; i++) {
		P[i]=[]; adjP[i]=[]; fitAdj[i+1]=[]

		adjP[i][0]=mol[i+1].atomicnumber
		for(j=1; j<NMAX; j++) adjP[i][j]=0

		for(j=0; j<n; j++) {
			k=adj[i+1][j+1]
			P[i][j]=k; fitAdj[i+1][j+1]=k
			adjP[i][mol[j+1].atomicnumber] += k
		}

		fitX[i+1]=mol[i+1].xini
		fitY[i+1]=mol[i+1].yini
		fitZ[i+1]=mol[i+1].zini
	}

	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			if(adjG[i][0]==1 || adjG[j][0]==1) { G[i][j]=0; G[j][i]=0}
			if(adjP[i][0]==1 || adjP[j][0]==1) { P[i][j]=0; P[j][i]=0}
		}
	}

	for(i=0; i<n; i++) {
		isMatch[i]=[]
		for(j=0; j<n; j++) {
			isMatch[i][j]=1
			for(k=0; k<NMAX; k++) {
				if(adjP[i][k]!=adjG[j][k]) {isMatch[i][j]=0; break}
			}
		}
	}

	var numAnch=0, anchG=[], anchP=[], txt=trim(norm($('idxFit').value))
	txt=txt.split('\n')

	for(k=0; k<txt.length; k++) {
		if(txt[k].indexOf('*')!=-1) {
			i=parseInt(txt[k].replace('*', ''))
			anchG[numAnch]=k
			anchP[numAnch]=i-1
			numAnch++
			for(j=0; j<n; j++) {
				isMatch[j][k]=0
				isMatch[i-1][j]=0
			}
			isMatch[i-1][k]=1
		}
	}

	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			if(isMatch[i][j]==1) {
				for(k=0; k<numAnch; k++) {
					if(bfs(P, anchP[k], i)!=bfs(G, anchG[k], j)) {isMatch[i][j]=0; break}
				}
			}
		}
	}

	t=getIsomorphicSubgraphs(G, P, 1);

	var idxFit=[]
	for(j=0; j<n; j++) {
		for(i=0; i<n; i++) if(t[0][i][j]==1) idxFit[j]=i+1
	}

	activeWin('molRef');
	mol=Mol(); n=mol[0].numatoms; adj=BondMatrix()
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) G[i][j]=adj[i+1][j+1]
	}

	activeWin('molFit');
	mol=Mol(); n=mol[0].numatoms; adj=BondMatrix()
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) P[i][j]=adj[i+1][j+1]
	}

	var used=[]
	for(i=0; i<n; i++) {
		if(adjP[i][0]==1) used[i]=false
		else used[i]=true
	}

	for(j=0; j<n; j++) {
		if(adjG[j][0]==1) {
			for(k=0; k<n; k++) if(G[j][k]) break
			for(i=0; i<n; i++) {
				if(!used[i] && adjP[i][0]==1 && P[i][idxFit[k]-1]) {
					idxFit[j]=i+1; used[i]=true; break
				}
			}
		}
	}
	$('idxFit').value=idxFit.join('\n')
}

function app() {
	idxFit='\n'+trim(norm($('idxFit').value))
	idxFit=idxFit.split('\n')

	activeWin('molFit');
	mol=Mol(); n=mol[0].numatoms; adj=BondMatrix()

	delMolecule();
	for(i=1; i<=n; i++) {
		idx=idxFit[i]
		addAtom(adjP[idx-1][0], fitX[idx], fitY[idx], fitZ[idx])
		for(j=1; j<i; j++) {
			jdx=idxFit[j]
			if(fitAdj[idx][jdx]) addBond(i, j)
		}
	}
	formula()
	centerMolecule();
	drawMolecule();
	InfoWin(idxFit.join(' '), 1)
}

function bfs(map, src, end){
	if(src==end) return 0

	let i, v, n=map.length, que=[], used=[], dist=[]

	i=n; while(i--) { dist[i]=0; used[i]=false }

	que.push(src)
	while(que.length) {
		v=que.shift()
		used[v]=true
		for(i=0; i<n; i++) {
			if( !used[i] && map[v][i]>0 ) {
				dist[i] = dist[v]+1
				if(i==end) return dist[i]
				used[i] = true
				que.push(i)
			}
		}
	}
	return null
}

function array2DCopy(A) {
    let X = [];
    for (let i = 0; i < A.length; i++) {
        X[i] = A[i].slice();
    }
    return X;
}

function mapPtoG(M) {
    return function (p) {
        const cols = M[0].length;
        for (let c = 0; c < cols; c++) {
            if (M[p][c] === 1) return c;
        }
    }
}

function isIso(M, G, P) {
    const rows = P.length;
    const morph = mapPtoG(M);

    for (let r1 = 0; r1 < rows; r1++) {
        for (let r2 = 0; r2 < rows; r2++) {
            // adjacent in P
            if (P[r1][r2] === 1) {
                // find mapped nodes in G
                let c1 = morph(r1);
                let c2 = morph(r2);
                // are they adjacent in G?
                if (G[c1][c2] !== 1) {
                    // no - not isomorphism
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 *
 * @param used_columns {number[]}
 * @param cur_row {number}
 * @param G {number[][]}
 * @param P {number[][]}
 * @param M {number[][]}
 * @param out {number[][][]}
 * @param num {number}
 * @param prune {boolean}
 */
function recurse(used_columns, cur_row, G, P, M, out, num, prune) {
    const cols = M[0].length;

    if (cur_row === M.length) {
        if (isIso(M, G, P)) out.push(array2DCopy(M))
    } else {
        let Mp = array2DCopy(M);

        // prune the proposed morphism to remove
        // mappings that are obviously not possible.
        if(prune){
            pruneOptions(Mp, P, G);
        }

        // for all unused columns c
        for (let c=0; c<cols; c++) {
            // only explore if the nodes are candidates for matching and the
            // column has not been set yet.
            if (used_columns[c] === 0 && M[cur_row][c] === 1) {

                //set column c in M' to 1 and other columns to 0
                for (let i=0; i<cols; i++) Mp[cur_row][i] = 0
                Mp[cur_row][c] = 1

                // mark c as used
                used_columns[c] = 1;

                // recurse, but only if they want to find more isomorphisms.
                if(num === null || out.length < num){
                    recurse(used_columns, cur_row+1, G, P, Mp, out, num);
                }

                // mark c as unused
                used_columns[c] = 0;
            }
        }
    }
}

/**
 *
 * @param M {number[][]} the proposed morphism between P and G
 * @param P {number[][]} the sub graph being matched
 * @param G {number[][]} the host graph
 */
function pruneOptions(M, P, G) {
    // M first dim (rows) are vertices of sub graph P
    // M second dim (cols) are vertices of host graph G
    console.log(M.length)
    for(let i = 0; i < M.length; i++){ // i - the vertex in P
        for(let j = 0; j < M.length; j++){ // j - the vertex in G
            if(M[i][j] === 1){ // for all M[i][j] === 1

                // for all neighbours x of vertex i in P
                for(let x = 0; x < P.length; x++){// x is a vertex in P that is adjacent to i
                    if(P[i][x] === 1){
                        // if there is no neighbour y of vertex j in G such
                        // that M[x][y] === 1, then set M[i][j] = 0

                        let hasNeighbourY = false;
                        for(let y = 0; y < G.length; y++){
                            if(G[j][y] === 1){
                                hasNeighbourY = true;
                                break;
                            }
                        }

                        if(!hasNeighbourY){
                            M[i][j] = 0;
                        }
                    }
                }
            }
        }
    }
}

/**
 * Finds isomorphisms (mappings) of a subgraph in a host/mother graph.
 *
 * The subgraph algorithm is based on: http://adriann.github.io/Ullman%20subgraph%20isomorphism.html
 *
 * This algorithm is exponential and will be slow for large inputs.
 *
 * @param G {number[][]} Adjacency matrix of the host/mother graph in which to search for a match.
 * @param P {number[][]} Adjacency matrix of subgraph to search for
 * @param maxNum {number} [null] the maximum number isomorphisms to find, may return fewer if fewer are matched.
 *
 * @returns {number[][][]} an array of morphism matrices (rows indices correspond to vertices of P, col indices correspond to vertices of G), null if error.
 */
function getIsomorphicSubgraphs(G, P, maxNum) {

    const G_size = G.length, P_size = P.length;

    // No match possible if |P| > |G|, not an error.
    // They don't want a match, not an error.
    if( (G_size < P_size) ||
		(maxNum !== null && maxNum !== undefined && maxNum <= 0) ) return [];

    let M = isMatch, used_col=[], results = [];

    i=G_size; while(i--) used_col[i]=0
    recurse(used_col, 0, G, P, M, results, maxNum, false);

    return results;
}

function saveFile(ext) {

	activeWin('molFit')
	var i, Lab, mol=Mol(), n=mol[0].numatoms

	var type='text', value
	if(ext=='pdb') { name='conf.pdb'
		value='TITLE Mol After Match\n'+'REMARK'+$('information').value+'\n'
		for(i=1; i<=n; i++) {
			Lab=trim(element(mol[i].atomicnumber,"symbol"))
			value += printf("%6s%5f%3s%16s%8.3f%8.3f%8.3f%24s\n",
				'HETATM', i, Lab, '0    ', mol[i].xini, mol[i].yini, mol[i].zini, Lab)
		}
		value += 'END\n'
	} else if(ext=='crd') { name='conf.inpcrd'
		value='MOL\n'+n+'\n'
		for(i=1; i<=n; i++) {
			value += printf("%12.7f%12.7f%12.7f\n", mol[i].xini, mol[i].yini, mol[i].zini)
		}
	}

	var blob;
	if (typeof window.Blob == "function") {
		blob = new Blob([value], {type: type});
	} else {
		var BlobBuilder =  window.BlobBuilder
						|| window.MozBlobBuilder
						|| window.WebKitBlobBuilder
						|| window.MSBlobBuilder;
		var bb = new BlobBuilder();
		bb.append(value);
		blob = bb.getBlob(type);
	}
	var URL = window.URL || window.webkitURL;
	var bloburl = URL.createObjectURL(blob);
	var anchor = document.createElement("a");
	if ('download' in anchor) {
		anchor.style.visibility = "hidden";
		anchor.href = bloburl;
		anchor.download = name;
		document.body.appendChild(anchor);
		var evt = document.createEvent("MouseEvents");
		evt.initEvent("click", true, true);
		anchor.dispatchEvent(evt);
		document.body.removeChild(anchor);
	} else if (navigator.msSaveBlob) {
		navigator.msSaveBlob(blob, name);
	} else {
		location.href = bloburl;
	}
}
//
// -------------------- periodic.js --------------------
//

//#   function element(Z,param)
//#
//#   Define public routine to store and return elemental information
//#
//#   Parameters:
//#      Z = Atomic number of element
//#      param:
//#         symbol   - Elemental symbol
//#         block    - "s", "p", "d", or "f" block
//#         valence  - # of valence electrons
//#         mass     - Atomic mass
//#         radius   - Covalent radius of atom
//#         EN       - Electronegativity of atom
//#         color    - color of atom (based on JMol)
//#         gradient - gradient color for atom
//#         label    - color for label
//#
//#   Return value:
//#      Variable containing requested information
//#
function element(Z,param) {
  if ( typeof element.elem == 'undefined' ) {
    element.elem = [];
    //
    //  Load parameters for all elements
    //
    //  Parameters:
    //    Atomic number
    //    Element symbol
    //    Block (s, p, d, or f)
    //    Number of valence electrons
    //    Atomic mass
    //    Covalent radii (in pm).  From www.webelements.com
    //    Jmol CPK colors (rgb).  From jmol.sourceforge.net/jscolors/
    //    Electronegativities from Wikipedia, citing WebElements as their source
    //
    element.elem[0]  = addElement(0, 'X', 's', 0, 0, 0, 0, 0, 0, 0);
    element.elem[1]  = addElement(1, 'H', 's', 1, 1.0079, 31, 255, 255, 255, 2.2);
    element.elem[2]  = addElement(2, 'He', 's', 2, 4.0026, 28, 217, 255, 255, 0);
    element.elem[3]  = addElement(3, 'Li', 's', 1, 6.941, 128, 204, 128, 255, 0.98);
    element.elem[4]  = addElement(4, 'Be', 's', 2, 9.0122, 96, 194, 255, 0, 1.57);
    element.elem[5]  = addElement(5, 'B', 'p', 3, 10.811, 84, 255, 181, 181, 2.04);
    element.elem[6]  = addElement(6, 'C', 'p', 4, 12.0107, 76, 144, 144, 144, 2.55);
    element.elem[7]  = addElement(7, 'N', 'p', 5, 14.0067, 71, 48, 80, 248, 3.04);
    element.elem[8]  = addElement(8, 'O', 'p', 6, 15.9994, 66, 255, 13, 13, 3.44);
    element.elem[9]  = addElement(9, 'F', 'p', 7, 18.9984, 57, 144, 224, 80, 3.98);
    element.elem[10] = addElement(10, 'Ne', 'p', 8, 20.1797, 58, 179, 227, 245, 0);
    element.elem[11] = addElement(11, 'Na', 's', 1, 22.9897, 166, 171, 92, 242, 0.93);
    element.elem[12] = addElement(12, 'Mg', 's', 2, 24.305, 141, 138, 255, 0, 1.31);
    element.elem[13] = addElement(13, 'Al', 'p', 3, 26.9815, 121, 191, 166, 166, 1.61);
    element.elem[14] = addElement(14, 'Si', 'p', 4, 28.0855, 111, 240, 200, 160, 1.9);
    element.elem[15] = addElement(15, 'P', 'p', 5, 30.9738, 107, 255, 128, 0, 2.19);
    element.elem[16] = addElement(16, 'S', 'p', 6, 32.065, 105, 255, 255, 48, 2.58);
    element.elem[17] = addElement(17, 'Cl', 'p', 7, 35.453, 102, 31, 240, 31, 3.16);
    element.elem[18] = addElement(18, 'Ar', 'p', 8, 39.948, 106, 128, 209, 227, 0);
    element.elem[19] = addElement(19, 'K', 's', 1, 39.0983, 203, 143, 64, 212, 0.82);
    element.elem[20] = addElement(20, 'Ca', 's', 2, 40.078, 176, 61, 255, 0, 1);
    element.elem[21] = addElement(21, 'Sc', 'd', 3, 44.9559, 170, 230, 230, 230, 1.36);
    element.elem[22] = addElement(22, 'Ti', 'd', 4, 47.867, 160, 191, 194, 199, 1.54);
    element.elem[23] = addElement(23, 'V', 'd', 5, 50.9415, 153, 166, 166, 171, 1.63);
    element.elem[24] = addElement(24, 'Cr', 'd', 6, 51.9961, 139, 138, 153, 199, 1.66);
    element.elem[25] = addElement(25, 'Mn', 'd', 7, 54.938, 139, 156, 122, 199, 1.55);
    element.elem[26] = addElement(26, 'Fe', 'd', 8, 55.845, 132, 224, 102, 51, 1.83);
    element.elem[27] = addElement(27, 'Co', 'd', 9, 58.9332, 126, 240, 144, 160, 1.88);
    element.elem[28] = addElement(28, 'Ni', 'd', 10, 58.6934, 124, 80, 208, 80, 1.91);
    element.elem[29] = addElement(29, 'Cu', 'd', 11, 63.546, 132, 200, 128, 51, 1.9);
    element.elem[30] = addElement(30, 'Zn', 'd', 2, 65.39, 122, 125, 128, 176, 1.65);
    element.elem[31] = addElement(31, 'Ga', 'p', 3, 69.723, 122, 194, 143, 143, 1.81);
    element.elem[32] = addElement(32, 'Ge', 'p', 4, 72.64, 120, 102, 143, 143, 2.01);
    element.elem[33] = addElement(33, 'As', 'p', 5, 74.9216, 119, 189, 128, 227, 2.18);
    element.elem[34] = addElement(34, 'Se', 'p', 6, 78.96, 120, 255, 161, 0, 2.55);
    element.elem[35] = addElement(35, 'Br', 'p', 7, 79.904, 120, 166, 41, 41, 2.96);
    element.elem[36] = addElement(36, 'Kr', 'p', 8, 83.8, 116, 92, 184, 209, 3);
    element.elem[37] = addElement(37, 'Rb', 's', 1, 85.4678, 220, 112, 46, 176, 0.82);
    element.elem[38] = addElement(38, 'Sr', 's', 2, 87.62, 195, 0, 255, 0, 0.95);
    element.elem[39] = addElement(39, 'Y', 'd', 3, 88.9059, 190, 148, 255, 255, 1.22);
    element.elem[40] = addElement(40, 'Zr', 'd', 4, 91.224, 175, 148, 224, 224, 1.33);
    element.elem[41] = addElement(41, 'Nb', 'd', 5, 92.9064, 164, 115, 194, 201, 1.6);
    element.elem[42] = addElement(42, 'Mo', 'd', 6, 95.94, 154, 84, 181, 181, 2.16);
    element.elem[43] = addElement(43, 'Tc', 'd', 7, 98, 147, 59, 158, 158, 1.9);
    element.elem[44] = addElement(44, 'Ru', 'd', 8, 101.07, 146, 36, 143, 143, 2.2);
    element.elem[45] = addElement(45, 'Rh', 'd', 9, 102.9055, 142, 10, 125, 140, 2.28);
    element.elem[46] = addElement(46, 'Pd', 'd', 10, 106.42, 139, 0, 105, 133, 2.2);
    element.elem[47] = addElement(47, 'Ag', 'd', 11, 107.8682, 145, 192, 192, 192, 1.93);
    element.elem[48] = addElement(48, 'Cd', 'd', 2, 112.411, 144, 255, 217, 143, 1.69);
    element.elem[49] = addElement(49, 'In', 'p', 3, 114.818, 142, 166, 117, 115, 1.78);
    element.elem[50] = addElement(50, 'Sn', 'p', 4, 118.71, 139, 102, 128, 128, 1.96);
    element.elem[51] = addElement(51, 'Sb', 'p', 5, 121.76, 139, 158, 99, 181, 2.05);
    element.elem[52] = addElement(52, 'Te', 'p', 6, 127.6, 138, 212, 122, 0, 2.1);
    element.elem[53] = addElement(53, 'I', 'p', 7, 126.9045, 139, 148, 0, 148, 2.66);
    element.elem[54] = addElement(54, 'Xe', 'p', 8, 131.293, 140, 66, 158, 176, 2.6);
    element.elem[55] = addElement(55, 'Cs', 's', 1, 132.9055, 244, 87, 23, 143, 0.79);
    element.elem[56] = addElement(56, 'Ba', 's', 2, 137.327, 215, 0, 201, 0, 0.89);
    element.elem[57] = addElement(57, 'La', 'd', 3, 138.9055, 207, 112, 212, 255, 1.1);
    element.elem[58] = addElement(58, 'Ce', 'f', 4, 140.116, 204, 255, 255, 199, 1.12);
    element.elem[59] = addElement(59, 'Pr', 'f', 5, 140.9077, 203, 217, 255, 199, 1.13);
    element.elem[60] = addElement(60, 'Nd', 'f', 6, 144.24, 201, 199, 255, 199, 1.14);
    element.elem[61] = addElement(61, 'Pm', 'f', 7, 145, 199, 163, 255, 199, 0);
    element.elem[62] = addElement(62, 'Sm', 'f', 8, 150.36, 198, 143, 255, 199, 1.17);
    element.elem[63] = addElement(63, 'Eu', 'f', 9, 151.964, 198, 97, 255, 199, 0);
    element.elem[64] = addElement(64, 'Gd', 'f', 10, 157.25, 196, 69, 255, 199, 1.2);
    element.elem[65] = addElement(65, 'Tb', 'f', 11, 158.9253, 194, 48, 255, 199, 0);
    element.elem[66] = addElement(66, 'Dy', 'f', 12, 162.5, 192, 31, 255, 199, 1.22);
    element.elem[67] = addElement(67, 'Ho', 'f', 13, 164.9303, 192, 0, 255, 156, 1.23);
    element.elem[68] = addElement(68, 'Er', 'f', 14, 167.259, 189, 0, 230, 117, 1.24);
    element.elem[69] = addElement(69, 'Tm', 'f', 15, 168.9342, 190, 0, 212, 82, 1.25);
    element.elem[70] = addElement(70, 'Yb', 'f', 16, 173.04, 187, 0, 191, 56, 0);
    element.elem[71] = addElement(71, 'Lu', 'f', 17, 174.967, 187, 0, 171, 36, 1.27);
    element.elem[72] = addElement(72, 'Hf', 'd', 4, 178.49, 175, 77, 194, 255, 1.3);
    element.elem[73] = addElement(73, 'Ta', 'd', 5, 180.9479, 170, 77, 166, 255, 1.5);
    element.elem[74] = addElement(74, 'W', 'd', 6, 183.84, 162, 33, 148, 214, 2.36);
    element.elem[75] = addElement(75, 'Re', 'd', 7, 186.207, 151, 38, 125, 171, 1.9);
    element.elem[76] = addElement(76, 'Os', 'd', 8, 190.23, 144, 38, 102, 150, 2.2);
    element.elem[77] = addElement(77, 'Ir', 'd', 9, 192.217, 141, 23, 84, 135, 2.2);
    element.elem[78] = addElement(78, 'Pt', 'd', 10, 195.078, 136, 208, 208, 224, 2.28);
    element.elem[79] = addElement(79, 'Au', 'd', 11, 196.9665, 136, 255, 209, 35, 2.54);
    element.elem[80] = addElement(80, 'Hg', 'd', 2, 200.59, 132, 184, 184, 208, 2);
    element.elem[81] = addElement(81, 'Tl', 'p', 3, 204.3833, 145, 166, 84, 77, 1.62);
    element.elem[82] = addElement(82, 'Pb', 'p', 4, 207.2, 146, 87, 89, 97, 2.33);
    element.elem[83] = addElement(83, 'Bi', 'p', 5, 208.9804, 148, 158, 79, 181, 2.02);
    element.elem[84] = addElement(84, 'Po', 'p', 6, 209, 140, 171, 92, 0, 2);
    element.elem[85] = addElement(85, 'At', 'p', 7, 210, 150, 117, 79, 69, 2.2);
    element.elem[86] = addElement(86, 'Rn', 'p', 8, 222, 150, 66, 130, 150, 0);
    element.elem[87] = addElement(87, 'Fr', 's', 1, 223, 260, 66, 0, 102, 0.7);
    element.elem[88] = addElement(88, 'Ra', 's', 2, 226, 221, 0, 125, 0, 0.9);
    element.elem[89] = addElement(89, 'Ac', 'd', 3, 227, 215, 112, 171, 250, 1.1);
    element.elem[90] = addElement(90, 'Th', 'f', 4, 232.0381, 206, 0, 186, 255, 1.3);
    element.elem[91] = addElement(91, 'Pa', 'f', 5, 231.0359, 200, 0, 161, 255, 1.5);
    element.elem[92] = addElement(92, 'U', 'f', 6, 238.0289, 196, 0, 143, 255, 1.38);
    element.elem[93] = addElement(93, 'Np', 'f', 7, 237, 190, 0, 128, 255, 1.36);
    element.elem[94] = addElement(94, 'Pu', 'f', 8, 244, 187, 0, 107, 255, 1.28);
    element.elem[95] = addElement(95, 'Am', 'f', 9, 243, 180, 84, 92, 242, 1.3);
    element.elem[96] = addElement(96, 'Cm', 'f', 10, 247, 169, 120, 92, 227, 1.3);
    element.elem[97] = addElement(97, 'Bk', 'f', 11, 247, 168, 138, 79, 227, 1.3);
    element.elem[98] = addElement(98, 'Cf', 'f', 12, 251, 168, 161, 54, 212, 1.3);
    element.elem[99] = addElement(99, 'Es', 'f', 13, 252, 165, 179, 31, 212, 1.3);
    element.elem[100] = addElement(100, 'Fm', 'f', 14, 257, 167, 179, 31, 186, 1.3);
    element.elem[101] = addElement(101, 'Md', 'f', 15, 258, 173, 179, 13, 166, 1.3);
    element.elem[102] = addElement(102, 'No', 'f', 16, 259, 176, 189, 13, 135, 1.3);
    element.elem[103] = addElement(103, 'Lr', 'f', 17, 262, 161, 199, 0, 102, 0);
    element.elem[104] = addElement(104, 'Rf', 'd', 4, 267, 157, 204, 0, 89, 0);
    element.elem[105] = addElement(105, 'Db', 'd', 5, 268, 149, 209, 0, 79, 0);
    element.elem[106] = addElement(106, 'Sg', 'd', 6, 269, 143, 217, 0, 69, 0);
    element.elem[107] = addElement(107, 'Bh', 'd', 7, 270, 141, 224, 0, 56, 0);
    element.elem[108] = addElement(108, 'Hs', 'd', 8, 269, 134, 230, 0, 46, 0);
    element.elem[109] = addElement(109, 'Mt', 'd', 9, 278, 129, 235, 0, 38, 0);
    element.elem[110] = addElement(110, 'Ds', 'd', 10, 281, 0, 0, 0, 28, 0);
    element.elem[111] = addElement(111, 'Rg', 'd', 11, 281, 0, 0, 0, 28, 0);
    element.elem[112] = addElement(112, 'Cn', 'd', 12, 285, 0, 0, 0, 28, 0);
    element.elem[113] = addElement(113, 'Uut', 'p', 3, 286, 0, 0, 0, 28, 0);
    element.elem[114] = addElement(114, 'Fl', 'p', 4, 289, 0, 0, 0, 28, 0);
    element.elem[115] = addElement(115, 'Uup', 'p', 5, 288, 0, 0, 0, 28, 0);
    element.elem[116] = addElement(116, 'Lv', 'p', 6, 293, 0, 0, 0, 28, 0);
    element.elem[117] = addElement(117, 'Uus', 'p', 7, 294, 0, 0, 0, 28, 0);
    element.elem[118] = addElement(118, 'Uuo', 'p', 8, 294, 0, 0, 0, 28, 0);
    }

  // Make sure element requested is defined.
  if ( (Z < 0)  ||  (Z > element.elem.length) ) {
    alert("Element "+Z+" is not defined.");
    return 0;
    }

  // Determine what information is requested and return value
  switch (param) {
    case "symbol":
        return element.elem[Z].symbol;
        break;
    case "block":
        return element.elem[Z].block;
        break;
    case "valence":
        return element.elem[Z].valence;
        break;
    case "mass":
        return element.elem[Z].mass;
        break;
    case "radius":
        return element.elem[Z].radius;
        break;
    case "EN":
        return element.elem[Z].EN;
        break;
    case "color":
        return element.elem[Z].color;
        break;
    case "gradient":
        return element.elem[Z].gradient;
        break;
    case "label":
        return element.elem[Z].label;
        break;
    case "max":
        return element.elem.length;
        break;
    default:
        alert("Param "+param+" is not defined.");
        return 0;
        break;
    }
  // End of element routine
  }

///
///   Define internal structure for elemental information
///
function elementObject() {
  this.symbol = "?";
  this.block = "?";
  this.valence = 0;
  this.mass = 0.0;
  this.radius = 0.0;
  this.EN = 0.0;
  this.color   = "rgb(88,88,88)";
  this.gradient = "rgb(88,88,88)";
  this.label = "rgb(255,255,255)";
  }

///
///   Internal routine to add information for elements
///
function addElement(Z, Symbol, block, valence, mass, radius, red, green, blue, EN)  {

  // Declare local variables
  var atomR, sum;
  var element = new elementObject();

  // Add elemental information
  element.symbol = Symbol;
  element.block = block;
  element.valence = valence;
  element.mass = mass;
  (radius <= 0) ? atomR = 0.1 : atomR = radius / 100.0;
  element.radius = atomR;
  element.color = "rgb("+red+","+green+","+blue+")";
  element.gradient = "rgb("+Math.floor(red/2)+","+Math.floor(green/2)+","+Math.floor(blue/2)+")";
  sum = red + green + blue;
  element.label = (sum > 384) ? "rgb(0,0,0)" : "rgb(255,255,255)";
  if(Z==6) element.label='rgb(255,255,255)'
  element.EN = EN;

  // End of addElement routine
  return element;
  }

//#   function drawPeriodic()
//#
//#   Writes html code to a division named "ptable".  The periodic table is written as a table, and each element
//#   is linked to to the pickElem() routine.
//#
//
// -------------------- molecule.js --------------------
//

//#   function Mol(value)
//#
//#   Define public routine to store and return molecule arrays
//#
//#   Parameter:
//#     value = If value given, then set as active molecule
//#
//#   Return value:
//#     Pointer to current array for molecular coordinates, ex.
//#        var molecule = Mol();
//#
//#   Molecular values stored in 'zeroth' element of array
//#     molecule[0].molIndex    - ID number for this molecule
//#     molecule[0].numatoms    - # of atoms in this molecule
//#     molecule[0].AtomScale   - Scale factor to control size of molecule displayed
//#     molecule[0].showlabels  - 0=don't show elemental symbols, 1=show symbols
//#     molecule[0].showcharges - 0=don't show charges, 1=show charges
//#     molecule[0].gradients   - 0=don't shade atoms, 1=use shading
//#     molecule[0].formula     - Simple text string containing molecular formula
//#     molecule[0].weight      - Molecular weight
//#     molecule[0].charge      - Charge of molecule
//#
//#   Atomic values in remaining elements of array
//#     molecule[i].atomicnumber - Atomic number (Z) for atom
//#     molecule[i].x            - x coordinate of atom
//#     molecule[i].y            - y coordinate of atom
//#     molecule[i].z            - z coordinate of atom
//#     molecule[i].charge       - Estimate of atomic charge
//#     molecule[i].highlite     - 0=don't highlight atom, 1=highlight atom
//#     molecule[i].hide         - 0=display atom, 1=do not display atom
//#
function Mol(value) {
  // Declare local variables
  var id;
  var param = value || 0;

  // Determine which molecule to use
  if(typeof Mol.molecule === "undefined")
    Mol.molecule = 1;
  if ( param > 0 ) {
    Mol.molecule = param;
    BondMatrix(param);
    }
  id = Mol.molecule;

  // Initialize molecular array if necessary
  if ( typeof Mol.moly === "undefined" )
    Mol.moly = [];
  if ( typeof Mol.moly[id] === "undefined" ) {
    Mol.moly[id] = [];
    Mol.moly[id][0] = new molObject();
    Mol.moly[id][0].molIndex = id;
    Mol.moly[id][0].numatoms = 0;
    Mol.moly[id][0].AtomScale = 1.0;
    Mol.moly[id][0].showlabels = 0;
    Mol.moly[id][0].showcharges = 0;
    Mol.moly[id][0].gradients = 1;
    Mol.moly[id][0].highlite = 0;
    Mol.moly[id][0].formula = "";
    Mol.moly[id][0].weight = 0.0;
    Mol.moly[id][0].charge = 999;
    }

  // End of Mol routine
  return Mol.moly[id];
  }

//#   function BondMatrix(value)
//#
//#   Define public routine to store and return Bond Matrix
//#
//#   Parameter:
//#      value = If value given, then set as active matrix
//#
//#   Return value:
//#      Pointer to current bond matrix, ex.
//#        var bonds = BondMatrix();
//#        bonds[i][j] then gives indication of bond between atoms "i" and "j"
//#
function BondMatrix(value) {
  // Declare local variables
  var id;
  var param = value || 0;

  // Determine which bond matrix to use
  if(typeof BondMatrix.id === "undefined")
    BondMatrix.id = 1;
  if ( param > 0 )
    BondMatrix.id = param;
  id = BondMatrix.id;

  // Initialize molecular array if necessary
  if ( typeof BondMatrix.bonds === "undefined" )
    BondMatrix.bonds = [];
  if ( typeof BondMatrix.bonds[id] === "undefined" ) {
    BondMatrix.bonds[id] = [];
    }

  // End of BondMatrix routine
  return BondMatrix.bonds[id];
  }

///
///  Private routine to store name for active window.
///
///  Parameter:
///    mywin = If name of window exists, then set as active window/molecule.
///            Otherwise, create new window with this name.
///
///  Return value:
///    Name of active Canvas/window
///
function activeWin(mywin) {

  // Declare local variables
  var i, active;
  var param = mywin || " ";

  // Initialize links
  if (typeof activeWin.mol === "undefined") {
    activeWin.mol = [];
    activeWin.mol[0] = "Dummy";
    }

  // Define defaults for active window and molecule
  if (typeof activeWin.win === "undefined") {
    var tag = document.getElementsByTagName("canvas");
    if ( tag[0] ) {
        activeWin.win = tag[0].id;
        activeWin.mol[1] = activeWin.win;
      } else {
        // alert("No canvas element found.");
        return;
      }
    }

  // User requested to set/view Window information
  if ( param.length > 1 )
    activeWin.win = param;

  // Find active molecule
  active=-1;
  for (i=0; i<activeWin.mol.length; i++)
    if ( activeWin.mol[i] == activeWin.win )
      active = i;

  // If requested window not found, add to list
  if ( active < 0 ) {
    active = activeWin.mol.length;
    activeWin.mol[active] = activeWin.win;
    }

  // Make correct molecule active
  Mol(active);

  // Define routines used to handle mouse events
  var canvas = $(activeWin.win);
  if ( canvas ) {
    canvas.onmousedown = MouseDown;
    document.onmouseup = MouseUp;
    canvas.onmousemove = MouseMove;
    canvas.addEventListener("mousewheel", MouseWheel, false);  // IE9,Chrome,Safari,Opera
    canvas.addEventListener("DOMMouseScroll", MouseWheel, false); // Firefox
    canvas.addEventListener("touchstart", MouseDown);
    canvas.addEventListener("touchend", MouseUp);
    canvas.addEventListener("touchmove", MouseMove);
    }

  // Finished with activeWin routine
  return activeWin.win;
  }

///
/// Define structure for each atom
///
function atomObject() {
  this.atomicnumber = 0;
  this.x = 0.0;
  this.y = 0.0;
  this.z = 0.0;
  this.xini = 0.0;
  this.yini = 0.0;
  this.zini = 0.0;
  this.charge = 0.0;
  this.highlite = 0;
  }

///
/// Define structure for molecular information
///
function molObject() {
  this.molIndex = 1;
  this.numatoms = 0;
  this.AtomScale = 1.0;
  this.showlabels = 1;
  this.showcharges = 0;
  this.gradients = 1;
  this.highlite = 0;
  this.hide = 0;
  this.formula = "";
  this.weight = 0.0;
  this.charge = 999;
  }

function readPDBfile(lines) {
  // Declare local variables
  var i, j, ncol, numatoms;
  var A, x, y, z;
  var txt;
  var current = new Array();
  var molecule = Mol();

  // Reset molecule object and contents
  delMolecule();

  // Read data
  InfoWin("",1);

	numatoms=0
	for(i=0; i<lines.length; i++) {
		txt=trim(lines[i].substr(0,6))
		current = lines[i].split(/\s+/)
		ncol=current.length
		if(txt=='ATOM' || txt=='HETATM') {
			numatoms++
			A = lookupSymbol(current[ncol-2]);
				if ( (A<1) || (A>118) ) {
				InfoWin("*** ERROR Reading PDB file. ***\n");
				return;
			}
			x = parseFloat(lines[i].substr(31,8));
			y = parseFloat(lines[i].substr(39,8));
			z = parseFloat(lines[i].substr(47,8));
			addAtom(A, x, y, z);
			//InfoWin(" "+current[0]+"("+A+")  "+x+", "+y+", "+z+"\n");
		} else if(txt=='CONECT') {
			for(j=2; j<ncol-1; j++) addBond(current[1],current[j])
		}
	}

  // Finished with readXYZfile routine
  formula();
  //InfoWin("Finished reading PDB input file.\n");
  }

///   function readXYZfile(lines)
///
///   Read xyz file - Format compatible with openbabel
///
///   lines:    Array holding content of .xyz file as lines
///
///   Data stored in molecule array and bond matrix
///
function readXYZfile(lines) {
  // Declare local variables
  var i, j, numatoms;
  var A, x, y, z;
  var dx, dy, dz, r;
  var bond;
  var mystr;
  var current = new Array();
  var Qtitle = new Array();
  var molecule = Mol();

  // Reset molecule object and contents
  delMolecule();

  // Read data
  InfoWin("",1);
  numatoms = parseInt(lines[0].replace(/\s+/,""));
  if ( isNaN(numatoms) ) {
    InfoWin("*** ERROR Reading XYZ file. ***\n");
    return;
    }
  Qtitle = lines[1];
  Qtitle = Qtitle.replace(/\s+$/,"");
  InfoWin(Qtitle+"\n");
  InfoWin("Looking for "+numatoms+" atoms.\n");

  // Read coordinates
  for (i=2; i < (numatoms+2); i++) {
    current = lines[i].replace(/\s+/g," ").replace(/^\s+/,"").split(' ');
    A = lookupSymbol(current[0]);
    if ( (A<1) || (A>118) ) {
      InfoWin("*** ERROR Reading XYZ file. ***\n");
      return;
      }
    x = parseFloat(current[1]);
    y = parseFloat(current[2]);
    z = parseFloat(current[3]);
    addAtom(A, x, y, z);
    InfoWin(" "+current[0]+"("+A+")  "+x+", "+y+", "+z+"\n");
    }

  // Create bonds
  for (i=1; i < numatoms; i++) {
    for (j=i+1; j <= numatoms; j++) {
      bond = 1.2 * (element(molecule[i].atomicnumber,"radius") + element(molecule[j].atomicnumber,"radius"));
      dx = molecule[i].x - molecule[j].x;
      dy = molecule[i].y - molecule[j].y;
      dz = molecule[i].z - molecule[j].z;
      r = Math.sqrt(dx*dx+dy*dy+dz*dz);
      if (r <= bond)
        addBond(i,j);
      }
    }

  // Finished with readXYZfile routine
  formula();
  InfoWin("Finished reading XYZ input file.\n");
  }

///   function readMOLfile(molfile,lines)
///
///   Read MDL .mol file (V2000 format)
///
///   molfile:  name of .mol file
///   lines:    Array holding content of .mol file
///
///   Data stored in molecule array and bond matrix
///
function readMOLfile(molfile,lines) {
  // Declare local variables
  var i, j;
  var x, y, z;
  var bond1, bond2;
  var A, symbol;
  var Qtitle = new Array();
  var molecule = Mol();

  // Write file information to information window
  var molname = escape(molfile.name);
  InfoWin("File Name = "+molname,1);
  InfoWin("\nFile Type = "+molfile.type);
  InfoWin("\nFile Size = "+molfile.size);

  // Reset molecule object and contents
  delMolecule();

  // Read data
  Qtitle = lines[0];
  Qtitle = Qtitle.replace(/\s+$/,"");
  var molAtoms = parseFloat(lines[3].substr(0,3));
  var molBonds = parseFloat(lines[3].substr(3,3));
  InfoWin("\n\n"+Qtitle);
  InfoWin("\n"+molAtoms+" atoms and "+molBonds+" bonds.");

  // Read coordinates
  j = molAtoms+4;
  for (i=4; i < j; i++) {
    x = parseFloat(lines[i].substr(0,10));
    y = parseFloat(lines[i].substr(10,10));
    z = parseFloat(lines[i].substr(20,10));
    symbol = lines[i].substr(31,3);
    symbol = symbol.replace(/\s+/,"");
    A = lookupSymbol(symbol);
    InfoWin("\n   "+symbol+"("+A+")  "+x+", "+y+", "+z);
    addAtom(A, x, y, z);
    }
  j = 4 + molAtoms + molBonds;
  for (i = 4+molAtoms; i < j; i++) {
    bond1 = parseFloat(lines[i].substr(0,3));
    bond2 = parseFloat(lines[i].substr(3,3));
    InfoWin("\n   Bond from "+bond1+" to "+bond2);
    addBond(bond1,bond2);
    }
  // Finished
  formula();
  }

///
///   Read GAMESS input file contain full set of x,y,z coordinates
///
///   lines:    Array holding content of .inp file
///
function readINPfile(lines) {
  // Declare local variables
  var i, j;
  var x, y, z;
  var bond1, bond2;
  var A, symbol;
  var numatoms;
  var Qtitle = new Array();
  var molecule = Mol();

  // Inform user that file being read
  InfoWin("--- Reading Gamess input file ---\n",1);

  // Reset molecule object and contents
  delMolecule();

  // Read title from file
  i = 0;
  lines[i].toLowerCase();
  while ( lines[i].toLowerCase().indexOf('$data') < 0 ) {
    i++;
    }
  i++;
  Qtitle = lines[i];
  Qtitle = Qtitle.replace(/\s+$/,"");
  InfoWin("\n\nTitle: [" + Qtitle + "]\n");

  // Advance to start of coordinates
  i++;
  if ( lines[i].toLowerCase().indexOf('c1') < 0 )
    i++;

  // Read coordinates
  i++;
  while ( lines[i].toLowerCase().indexOf('$end') < 0 ) {
    InfoWin("   Line = [" + lines[i] + "]\n");
    current = lines[i].replace(/\s+/g," ").replace(/^\s+/,"").split(' ');
    A = parseInt(current[1]);
    x = parseFloat(current[2]);
    y = parseFloat(current[3]);
    z = parseFloat(current[4]);
    addAtom(A, x, y, z);
    i++;
    }

  // Calculate bonds
  for (i=1; i < molecule[0].numatoms; i++) {
    ra = element(molecule[i].atomicnumber,"radius");
    for (j=i+1; j <= molecule[0].numatoms; j++) {
      rb = element(molecule[i].atomicnumber,"radius");
      r = 0.9*distance(molecule, i, j);
      if ( r <= (ra+rb) )
        addBond(i, j);
      }
    }

  // Finished with readINPfile routine
  InfoWin("Finished reading input file.\n");
  formula();
  return;
  }

//#   function lookupSymbol(symbol)
//#
//#   Given elemental symbol, lookup and return atomic number (Z)
//#
function lookupSymbol(symbol) {
  var i;

  for (i = 1; i < element(1,"max"); i++) {
    if ( element(i,"symbol") == symbol )
      return i;
    }
  return 0;
  }

///
///   Routine to remove and initialize molecule object
///
function delMolecule() {
  var molecule = Mol();

  molecule[0] = new molObject();
  molecule[0].molIndex = 1;
  molecule[0].numatoms = 0;
  molecule[0].AtomScale = 1.0;
  molecule[0].showlabels = 0;
  molecule[0].showcharges = 0;
  buttonColor("ChargeButton",0);
  }

///
///  Routine to save/restore molecular coordinates, bonds, ...
///
///  Parameter: mode
///      pop  - restore last saved set of coordinates
///      save - store current information
///      redo - 'Undo' last pop
///
function Undo(mode) {

  // Declare local variables
  var i, pos;
  var maxBonds;
  var MAXSAVE=10;
  var molecule = Mol();
  var bonds = BondMatrix();

  // Initialization of variables
  if ( typeof Undo.num == 'undefined' ) {
    Undo.num = 0;
  	Undo.mol = [];
  	Undo.bonds = [];
  	for (i=0; i < MAXSAVE; i++) {
  	  Undo.mol[i] = [];
  	  Undo.mol[i][0] = 0;
  	  Undo.bonds[i] = [];
  	  }
    }

  // Save current geometry information
  if ( mode == "save")  {
  	pos = Undo.num;
    // InfoWin(" --- Saving geometry in slot " + pos + " ---\n");
  	Undo.mol[pos] = null;
  	Undo.mol[pos] = [];
    Undo.mol[pos][0] = molecule[0].numatoms;
    for (i=1; i <= molecule[0].numatoms; i++) {
      Undo.mol[pos][i] = new atomObject();
      Undo.mol[pos][i].atomicnumber = molecule[i].atomicnumber;
      Undo.mol[pos][i].x = molecule[i].x;
      Undo.mol[pos][i].y = molecule[i].y;
      Undo.mol[pos][i].z = molecule[i].z;
      Undo.mol[pos][i].bonds = [];
      Undo.mol[pos][i].charge = molecule[i].charge;
      Undo.mol[pos][i].highlite = molecule[i].highlite;
      }
  	Undo.bonds[pos] = null;
  	Undo.bonds[pos] = [];
    for (i=1; i <= molecule[0].numatoms; i++) {
      Undo.bonds[pos][i] = [];
      for (j=0; j <= molecule[0].numatoms; j++) {
        Undo.bonds[pos][i][j] = bonds[i][j];
        }
      }

    // alert("Coordinate set "+Undo.num+" contains "+Undo.mol[pos][0]+" atoms.");
    Undo.num = ( (Undo.num+1) % MAXSAVE );
    return;
    }

  // Set position of 'last saved' geometry information
  if ( mode == "pop")  {
    Undo('save');
    pos = ( (Undo.num+MAXSAVE-2) % MAXSAVE );
  	Undo.num = pos;
    // InfoWin(" --- Attempting to restore from slot " + pos + " ---\n");
  	if ( Undo.mol[pos][0] < 1 ) {
      alert("ERROR: No saved coordinates to restore.");
      return;
  	  }
  	}

  // Attempt to reverse last 'undo' (pop)
  if ( mode == "redo")  {
    pos = ( (Undo.num+1) % MAXSAVE );
  	Undo.num = pos;
    // InfoWin(" --- Attempting to go back to slot " + pos + " ---\n");
  	if ( Undo.mol[pos][0] < 1 ) {
      alert("ERROR: No saved coordinates to restore.");
      return;
  	  }
  	}

  if ( (mode == "redo") || (mode == "pop") )  {
    molecule[0].numatoms = Undo.mol[pos][0];
    for (i=1; i <= Undo.mol[pos][0]; i++) {
      molecule[i].atomicnumber = Undo.mol[pos][i].atomicnumber;
      molecule[i].x = Undo.mol[pos][i].x;
      molecule[i].y = Undo.mol[pos][i].y;
      molecule[i].z = Undo.mol[pos][i].z;
      molecule[i].charge = Undo.mol[pos][i].charge;
      molecule[i].highlite = Undo.mol[pos][i].highlite;
      }
    for (i=1; i <= molecule[0].numatoms; i++) {
      for (j=0; j <= molecule[0].numatoms; j++) {
        bonds[i][j] = Undo.bonds[pos][i][j];
        }
      }
    centerMolecule();
    drawMolecule();
    return;
    }

  // End of Undo routine
  }

//#   function addAtom(AtomicNum, x, y, z)
//#
//#   Routine to add atomic information to molecule array.
//#
//#   Parameters are the atomic # and (x,y,z) coordinates of the atom.
//#
function addAtom(AtomicNum, x, y, z) {

  // Declare local variables
  var i;
  var numatoms;
  var molecule = Mol();
  var bonds = BondMatrix();

  // Create/fill molecular array
  molecule[0].numatoms++;
  numatoms = molecule[0].numatoms;
  molecule[0].showcharges = 0;
  buttonColor("ChargeButton",0);
  molecule[numatoms] = new atomObject();
  molecule[numatoms].atomicnumber = AtomicNum;
  molecule[numatoms].x = x;
  molecule[numatoms].y = y;
  molecule[numatoms].z = z;
  molecule[numatoms].xini = x;
  molecule[numatoms].yini = y;
  molecule[numatoms].zini = z;
  molecule[numatoms].charge = 0.0;
  molecule[numatoms].highlite = 0;

  // Initialize bond matrix
  if (typeof bonds[numatoms] === "undefined") {
    bonds[numatoms] = [];
    }
  bonds[numatoms][0] = 0;
  for (i=1; i <= numatoms; i++) {
    bonds[i][numatoms] = 0;
    bonds[numatoms][i] = 0;
    }

  // End of addAtom routine
  }

//#   function addBond(atom1, atom2)
//#
//#   Routine to add bonding information to bonds array
//#   (Note that atom numbers start at 1)
//#
function addBond(atom1, atom2) {

  // Define local variables
  var i, j;
  var newbond;
  var molecule = Mol();
  var numatoms = molecule[0].numatoms;
  var bonds = BondMatrix();

  // Add bonds to matrix
  bonds[atom1][0]++;
  bonds[atom1][atom2] = 1;
  bonds[atom2][0]++;
  bonds[atom2][atom1] = 1;

  // End of addBond routine
  }

//#   function delAtom(atomNum)
//#
//#   Routine to remove selected atom, directly-bonded Hydrogens, and
//#   cleanup molecule object (remove blank slots) and update bond arrays.
//#
//#   atomNum is the position of the atom to delete in the molecule array
//#
function delAtom(atomNum) {

  // Declare local variables
  var i, j, k, Hi=0;
  var bondedAtom, bi;
  var numatoms;
  var newNum, maxBonds, maxBondsB, maxBondsI, maxBondsJ;
  var HList = new Array();
  var molecule = Mol();
  var bonds = BondMatrix();

  // See if there are any attached hydrogen atoms
  molecule[0].showcharges = 0;
  buttonColor("ChargeButton",0);
  newNum = atomNum;
  for (i=1; i < atomNum; i++) {
    if ( (bonds[atomNum][i] > 0) && (molecule[i].atomicnumber == 1) ) {
      HList[Hi++] = i;
      newNum--;
      }
    }
  for (i=atomNum+1; i <= molecule[0].numatoms; i++) {
    if ( (bonds[atomNum][i] > 0) && (molecule[i].atomicnumber == 1) ) {
      HList[Hi++] = i;
      }
    }

  // Remove attached hydrogen atoms
  HList.sort(function(a,b){
    return(b-a);
    });
  for (i=0; i < Hi; i++) {
  	delAtom(HList[i]);
    }
  atomNum = newNum;

  // Remove atom and shift 'higher' atoms in molecule object to fill gaps
  for (i=atomNum; i < molecule[0].numatoms; i++) {
  	j=i+1;
    molecule[i].atomicnumber = molecule[j].atomicnumber;
    molecule[i].x = molecule[j].x;
    molecule[i].y = molecule[j].y;
    molecule[i].z = molecule[j].z;
    molecule[i].charge = molecule[j].charge;
    molecule[i].highlite = molecule[j].highlite;
    }
  molecule[0].numatoms--;

  // Remove bonds to atom being deleted
  numatoms = molecule[0].numatoms;
  for (i=1; i < atomNum; i++) {
    if ( bonds[i][atomNum] > 0 )
      bonds[i][0]--;
    for (j=atomNum; j <= numatoms; j++)
      bonds[i][j] = bonds[i][j+1];
    bonds[i][numatoms+1] = 0;
    }
  for (i=atomNum; i <= numatoms; i++) {
    if ( bonds[i+1][atomNum] > 0 )
      bonds[i+1][0]--;
    for (j=0; j < atomNum ; j++)
      bonds[i][j] = bonds[i+1][j];
    for (j=atomNum; j <= numatoms; j++)
      bonds[i][j] = bonds[i+1][j+1];
    bonds[i][numatoms+1] = 0;
    }
  for (i=0; i<=numatoms+1; i++)
    bonds[numatoms+1][i] = 0;

  // Finished with delAtom routine
  drawMolecule();
  }

//#   function delBond(atom1, atom2)
//#
//#   Routine to remove bond from bonds array.
//#
function delBond(atom1, atom2) {
  var i, maxBonds;
  var bonds = BondMatrix();

  // Remove bond on atom1
  if ( bonds[atom1][atom2] > 0 ) {
    bonds[atom1][atom2] = 0;
    bonds[atom1][0]--;
    }
  // Remove bond on atom2
  if ( bonds[atom2][atom1] > 0 ) {
    bonds[atom2][atom1] = 0;
    bonds[atom2][0]--;
    }

  // Finished with delBond routine
  }

//#   function hideH()
//#
//#   Do not display any hydrogen atoms in molecule
//#
function hideH() {
  var i;
  var molecule = Mol();

  // Loop over all atoms, setting hide flag to 0 for all atoms
  for (i=1; i <= molecule[0].numatoms; i++) {
    if ( molecule[i].atomicnumber == 1 ) {
      molecule[i].hide = 1;
      }
    }

  // Finished with hideH routine
  drawMolecule();
  }

//#   function showAll()
//#
//#   Show all atoms in molecule
//#
function showAll() {
  var i;
  var molecule = Mol();

  // Loop over all atoms, setting hide flag to 0 for all atoms
  for (i=1; i < molecule[0].numatoms; i++) {
    molecule[i].hide = 0;
    }

  // Finished with showAll routine
  }

//#   function centerMolecule()
//#
//#   Find center of molecule and move coordinates to center.
//#
function centerMolecule() {

  // Define local variables
  var i, A, numatoms;
  var molSize = 0.0;
  var xc, yc, zc;
  var xyz;
  var activeCanvas = activeWin("");
  var molecule = Mol();

  // Get height of Canvas
  var canvas = $(activeCanvas);
  var height = canvas.height;

  // Find center of molecule
  xc = 0.0;
  yc = 0.0;
  zc = 0.0;
  numatoms = molecule[0].numatoms;
  for (i=1; i <= numatoms; i++) {
    xc += molecule[i].x;
    yc += molecule[i].y;
    zc += molecule[i].z;
    }
  xc /= numatoms;
  yc /= numatoms;
  zc /= numatoms;

  // Translate coordinates to center
  molSize = 0.0;
  for (i=1; i <= numatoms; i++) {
    A = molecule[i].atomicnumber;
    molecule[i].x -= xc;
    molecule[i].y -= yc;
    molecule[i].z -= zc;
    xyz = Math.abs(molecule[i].x);
    xyz = Math.max(xyz,Math.abs(molecule[i].y));
    xyz = Math.max(xyz,Math.abs(molecule[i].z));
    xyz += element(A,"radius");
    molSize = Math.max(molSize,xyz);
    }
  molecule[0].AtomScale = 0.5*height/molSize;

  // Finished with centerMolecule routine
  }

//#   function formula()
//#
//#   Display molecular formula to division with an id = "formula".
//#   Note that the contents of this division will be overwritten.
//#
function formula() {
  var i, j, A;
  var found, num;
  var molweight;
  var molformula = new Array();
  var formulaStr="Formula = ";
  var molecule = Mol();

  // If no molecule loaded, return
  if (molecule[0].numatoms < 1)
    return;

  // Look for carbon and place in first slot
  num = 0;
  for (i=1; i<=molecule[0].numatoms; i++)
    if ( molecule[i].atomicnumber == 6 ) {
      molformula[num] = [];
      molformula[num][0] = 6;
      molformula[num][1] = 0;
      num++;
      i = molecule[0].numatoms+1;
      }

  // Look for hydrogen in next slot
  for (i=1; i<=molecule[0].numatoms; i++)
    if ( molecule[i].atomicnumber == 1 ) {
      molformula[num] = [];
      molformula[num][0] = 1;
      molformula[num][1] = 0;
      num++;
      i = molecule[0].numatoms+1;
      }

  // Loop over all atoms in molecule
  for (i=1; i<=molecule[0].numatoms; i++) {
    A = molecule[i].atomicnumber;
    found = 0;
    for (j=0; j < num; j++) {
      if (molformula[j][0] == A) {
        molformula[j][1]++;
        found = 1;
        j = num;
        }
      }
    if (found == 0) {
      molformula[num] = new Array(2);
      molformula[num][0] = molecule[i].atomicnumber;
      molformula[num][1] = 1;
      num++;
      }
    }

  // Calculate molecular weight and output formula
  molweight = 0.0;
  molecule[0].formula = "";
  for (i=0; i<num; i++) {
    A = molformula[i][0];
    molweight += element(A,"mass") * molformula[i][1];
    formulaStr += element(A,"symbol");
    molecule[0].formula += element(A,"symbol");
    if (molformula[i][1] > 1) {
      formulaStr += "<sub>"+molformula[i][1]+"</sub>";
      molecule[0].formula += molformula[i][1] + " ";
      }
    }
  molecule[0].weight = molweight;

  // Write formula to screen
  if (activeWin()=='molRef') $('forRef').innerHTML = formulaStr;
  if (activeWin()=='molFit') $('forFit').innerHTML = formulaStr;
  if ( molweight > 0.0 )
    InfoWin(molecule[0].formula+" has a molecular weight of "+molweight.toFixed(2)+"\n");
  // End formula routine
	showLabels()
	var str=''
	for (i=1; i<=molecule[0].numatoms; i++) str += i+'\n'

  if (activeWin()=='molRef') $('idxRef').value = str;
  if (activeWin()=='molFit') $('idxFit').value = str;

  }

//#   function RotateMolecule(axis)
//#
//#   Routine to start/stop rotation of molecule.
//#
//#   Parameter:
//#     axis = rotation axis.  Allowed values are "x", "y", or "z"
//#            ("s" can be used to stop all rotations)
//#
//#   Buttons must be created in html with:
//#     id="rotateX"  (or rotateY or rotateZ)
//#

  function RotateMolecule(axis) {
    // Declare local variables
    var delay = 20;  // Delay between redrawing screen, in ms

    // Parameter validation
    axis = axis[0].toLowerCase();
    if ( (axis != "x") &&  (axis != "y") && (axis != "z") && (axis != "s") ) {
      return;
      }

    // If first call, initialize "active" variable
    if(typeof RotateMolecule.RX === "undefined") {
      RotateMolecule.RX = 0;
      RotateMolecule.RY = 0;
      RotateMolecule.RZ = 0;
      }

    // Start/Stop rotation
    switch (axis) {
      case 'x':
                if ( RotateMolecule.RX == 0 ) {
                    // Start rotation
                    RotateMolecule.RX = 1;
                    RotateX_ID = setInterval( "rotateX()", delay );
                    buttonColor("rotateX",1);
                  } else {
                    // Stop rotation
                    RotateMolecule.RX = 0;
                    clearInterval( RotateX_ID );
                    buttonColor("rotateX",0);
                  }
                break;
      case 'y':
                if ( RotateMolecule.RY == 0 ) {
                    // Start rotation
                    RotateMolecule.RY = 1;
                    RotateY_ID = setInterval( "rotateY()", delay );
                    buttonColor("rotateY",1);
                  } else {
                    // Stop rotation
                    RotateMolecule.RY = 0;
                    clearInterval( RotateY_ID );
                    buttonColor("rotateY",0);
                  }
                break;
      case 'z':
                if ( RotateMolecule.RZ == 0 ) {
                    // Start rotation
                    RotateMolecule.RZ = 1;
                    RotateZ_ID = setInterval( "rotateZ()", delay );
                    buttonColor("rotateZ",1);
                  } else {
                    // Stop rotation
                    RotateMolecule.RZ = 0;
                    clearInterval( RotateZ_ID );
                    buttonColor("rotateZ",0);
                  }
                break;
      case 's':
                if ( RotateMolecule.RX != 0 ) {
                  RotateMolecule.RX = 0;
                  clearInterval( RotateX_ID);
                  }
                if ( RotateMolecule.RY != 0 ) {
                  RotateMolecule.RY = 0;
                  clearInterval( RotateY_ID);
                  }
                if ( RotateMolecule.RZ != 0 ) {
                  RotateMolecule.RZ = 0;
                  clearInterval( RotateZ_ID);
                  }
                buttonColor("rotateX",0);
                buttonColor("rotateY",0);
                buttonColor("rotateZ",0);
                break;
      }

    // End of RotateMolecule function
    }

//
//   Private function to perform rotation of molecule about x-axis
//
  function rotateX() {
    // Declare local variables
    var i, x, y, z;
    var angle = 1.0;
    var cosA = Math.cos(-angle*Math.PI/180);
    var sinA = Math.sin(-angle*Math.PI/180);
    var molecule = Mol();

    // Draw molecule and display on canvas
    for (i=1; i <= molecule[0].numatoms; i++) {
      y = molecule[i].y;
      z = molecule[i].z;
      molecule[i].y =  cosA*y + sinA*z;
      molecule[i].z = -sinA*y + cosA*z;
      }
    drawMolecule();

    // End of rotateX function
    }

//
//   Private function to perform rotation of molecule about y-axis
//
  function rotateY() {
    // Declare local variables
    var i, x, y, z;
    var angle = 1.0;
    var cosA = Math.cos(-angle*Math.PI/180);
    var sinA = Math.sin(-angle*Math.PI/180);
    var molecule = Mol();

    // Draw molecule and display on canvas
    for (i=1; i <= molecule[0].numatoms; i++) {
      x = molecule[i].x;
      z = molecule[i].z;
      molecule[i].x =  cosA*x + sinA*z;
      molecule[i].z = -sinA*x + cosA*z;
      }
    drawMolecule();

    // End of rotateY function
    }

//
//   Private function to perform rotation of molecule about z-axis
//
  function rotateZ() {
    // Declare local variables
    var i, x, y, z;
    var angle = 1.0;
    var cosA = Math.cos(-angle*Math.PI/180);
    var sinA = Math.sin(-angle*Math.PI/180);
    var molecule = Mol();

    // Draw molecule and display on canvas
    for (i=1; i <= molecule[0].numatoms; i++) {
      x = molecule[i].x;
      y = molecule[i].y;
      molecule[i].x =  cosA*x + sinA*y;
      molecule[i].y = -sinA*x + cosA*y;
      }
    drawMolecule();

    // End of rotateZ function
    }
//
// -------------------- drawmolecule.js --------------------
//

//#   function drawMolecule()
//#
//#   Draw atoms and bonds.
//#
function drawMolecule() {

  // Define local variables
  var i, j, k;
  var x, y, r;
  var bonds = BondMatrix();
  var molecule = Mol();
  var A, numatoms=molecule[0].numatoms;
  var bondto, maxBonds;
  var x1, y1, x2, y2;
  var AtomSize = 0.50;
  var SimpleSwitch = 200;

  // Create connection to active Canvas
  var mycanvas = activeWin("");
  var canvas = $(mycanvas);
  if ( ! canvas )
    return;
  var ctx = canvas.getContext('2d');
  var width = canvas.width;
  var height = canvas.height;
  var centerx = canvas.width/2;
  var centery = canvas.height/2;

  // Clear screen before drawing molecule
  clear(ctx, width, height);

  // Define properties of lines
  ctx.strokeStyle = "rgb(0, 0, 255)";
  ctx.lineWidth = 3.0;

  // Depth sort: Draw from back to front
  if ( typeof drawMolecule.deep == 'undefined' ) {
    drawMolecule.deep = [];
    }
  if ( drawMolecule.deep.length != numatoms ) {
    j=0;
    drawMolecule.deep = [];
    for (i=drawMolecule.deep.length; i<numatoms; i++) {
      j++;
      drawMolecule.deep[i] = [];
      drawMolecule.deep[i]["id"] = j;
      drawMolecule.deep[i]["sz"] = molecule[j].z;
      }
    }
  for (i=0; i<numatoms; i++) {
    j = drawMolecule.deep[i].id;
    drawMolecule.deep[i].z = molecule[j].z;
    }
  drawMolecule.deep.sort(function(a,b){
    return(b.z - a.z);
    });

  // If too many atoms, do not use gradients
  if ( numatoms > SimpleSwitch )
    molecule[0].gradients = 0;

  // Loop over all atoms in molecule
  for (i = 0; i < numatoms; i++) {
    atom = drawMolecule.deep[i].id;
    if ( ! molecule[atom].hide ) {
          drawAtom(ctx, atom, AtomSize, centerx, centery);
      if ( molecule[0].showcharges != 0 )
        DrawChargeCloud(ctx, atom, AtomSize*1.5, centerx, centery);
      if ( molecule[atom].highlite != 0 )
        DrawAtomHilite(ctx, atom, AtomSize*1.5, centerx, centery);
      // Loop over all bonds on current atom
      for (j = 1; j <= numatoms; j++) {
        if ( (bonds[atom][j] > 0) && (j != atom) && !molecule[j].hide ) {
          for (k=i+1; k < numatoms; k++) {
            if (drawMolecule.deep[k].id == j) {
              drawBond(ctx, atom, j, AtomSize, centerx, centery, molecule[0].gradients);
              }
            }
          }
        }
      }
    }

  // Finished with drawMolecule routine
  }

//#   function showLabels()
//#
//#   Routine to toggle display of atomic labels.
//#
function showLabels() {

  // Local variables
  var molecule = Mol();

  if ( molecule[0].showlabels == 0 ) {
    buttonColor("LabelButton",1);
    molecule[0].showlabels = 1;
    drawMolecule();
    return;
    }
  buttonColor("LabelButton",0);
  molecule[0].showlabels = 0;
  drawMolecule();
  return;
  }

///
///   Routine to draw a single atom
///
function drawAtom(ctx, atomNum, AtomSize, centerx, centery) {
  // Define local variables
  var A;
  var x, y, r, width;
  var xoff, yoff, roff;
  var molecule = Mol();
  var activeCanvas = activeWin("");
  // Get parameters for atom
  A = molecule[atomNum].atomicnumber;
  x = molecule[0].AtomScale*molecule[atomNum].x + centerx;
  y = molecule[0].AtomScale*molecule[atomNum].y + centery;
  r = molecule[0].AtomScale*AtomSize*element(A,"radius");
  xoff = x-(0.20*r);
  yoff = y-(0.20*r);
  roff = .30*r;
  x = Math.floor(x+0.5);
  y = Math.floor(y+0.5);

  if(x<0 || x>2*centerx || y<0 || y>2*centery) return

	if(drawMode<=1) {
		width=1
		r = Math.floor(r+0.5);
		AtomColor = element(A,"color")
		if(drawMode==1) {
			width=0.01
			xoff = Math.floor(xoff+0.5);
			yoff = Math.floor(yoff+0.5);
			AtomColor = ctx.createRadialGradient(xoff-0.5, yoff-0.5, roff, x, y, r);
			AtomColor.addColorStop(0, element(A,"color"));
			AtomColor.addColorStop(1, element(A,"gradient"));
		}

		// Plot atom as a gradient shaded circle
		ctx.beginPath();
		ctx.arc(x, y, r, 0, 2*Math.PI, false);
		ctx.fillStyle = AtomColor;
		ctx.fill();
		ctx.lineWidth = width;
		ctx.strokeStyle = "black";
		ctx.stroke();
		ctx.closePath();
	} else if(drawMode==2) ctx.drawImage(img[A], x-r, y-r, r*2, r*2);

  // Draw label for atom
	if ( molecule[0].showlabels != 0 ) atomLabel(ctx,A,atomNum,x,y);
}
  // Plot atom as a solid shaded circle
///
///   Routine to draw a single atom without gradient
///

///
///   Routine to draw a single atom with a highlight
///
function DrawAtomHilite(ctx, atomNum, AtomSize, centerx, centery) {

  // Define local variables
  var A;
  var x, y, r;
  var colR, colB, colG, deltaRGB;
  var molecule = Mol();
  var activeCanvas = activeWin("");

  // Get parameters for atom
  A = molecule[atomNum].atomicnumber;
  x = molecule[0].AtomScale*molecule[atomNum].x + centerx;
  y = molecule[0].AtomScale*molecule[atomNum].y + centery;
  r = molecule[0].AtomScale*AtomSize*element(A,"radius");

  // Set highlight color to yellow
  AtomColor = "rgba(255,255,0,0.5)";

  // Plot charge cloud as a semi-transparent circle
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI, false);
  ctx.fillStyle = AtomColor;
  ctx.fill();
  ctx.lineWidth = 0.01;
  ctx.strokeStyle = "black";
  ctx.stroke();
  ctx.closePath();

  // End DrawAtomHilite routine
  }

///
///   Helper routine to draw label for atom
///   Parameters:
///     Atomic number, atom #, (x,y) coordinates
///
function atomLabel(ctx,A,i,ax,ay) {
  var label = element(A,"symbol")+i;
  var fontsize = 14;
  var x = parseFloat(ax);
  var y = parseFloat(ay) + fontsize/2;
  var molecule = Mol();

  ctx.lineWidth = 1.0;
  ctx.textAlign = "center";
  ctx.font = "normal " + fontsize + "px sans-serif";
  ctx.fillStyle = element(A,"label");
  ctx.beginPath();
  ctx.fillText(label,x,y);
  var w = 10.0 / molecule[0].AtomScale;
  ctx.fillRect(x,y,w*1.5,w);
  ctx.closePath();
  }

///
///   Routine to draw a single bond
///
function drawBond(ctx, atom1, atom2, AtomSize, centerx, centery, smallMolecule) {

  // Define local variables
  var i, j, k;
  var x1, y1, z1, x2, y2, z2;
  var dx, dy, dz;
  var xoff, yoff;
  var dist, r, width;
  var molecule = Mol();
  var activeCanvas = activeWin("");
  var AtomScale = molecule[0].AtomScale;
  var BondWidth = Math.floor(0.10*AtomScale + 0.5);
  var BondColor = "rgb(200,128,51)";

  // Get parameters for atoms
  r =  AtomScale * AtomSize * element(molecule[atom1].atomicnumber,"radius") - 2;
  x1 = AtomScale*molecule[atom1].x + centerx;
  y1 = AtomScale*molecule[atom1].y + centery;
  z1 = AtomScale*molecule[atom1].z;
  x2 = AtomScale*molecule[atom2].x + centerx;
  y2 = AtomScale*molecule[atom2].y + centery;
  z2 = AtomScale*molecule[atom2].z;
  x1 = Math.floor(x1+0.5);
  y1 = Math.floor(y1+0.5);
  x2 = Math.floor(x2+0.5);
  y2 = Math.floor(y2+0.5);

  // Calculate intersection of bond with sphere
  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  dist = Math.sqrt(dx*dx + dy*dy + dz*dz);

  xoff = x1 + r*dx/dist;
  yoff = y1 + r*dy/dist;
	xoff=Math.floor(xoff+0.5);
	yoff=Math.floor(yoff+0.5);

  // For large molecules, use simpler plot
	// Draw bond
	if(drawMode<=1) {
		ctx.beginPath();
		ctx.moveTo(xoff, yoff);
		ctx.lineTo(x2, y2);

		ctx.lineWidth = BondWidth;
		ctx.strokeStyle = BondColor;
		ctx.stroke();

		if(drawMode==1) {
			ctx.lineWidth = BondWidth + 2;
			ctx.strokeStyle = "rgb(1,0,0)";
			ctx.stroke();
			ctx.moveTo(xoff, yoff);
			ctx.lineTo(x2, y2);
			ctx.lineWidth = BondWidth;
			ctx.strokeStyle = BondColor;
			ctx.stroke();
		}
		ctx.closePath();
	} else if(drawMode==2) {
		r=Math.sqrt((y2-yoff)*(y2-yoff)+(x2-xoff)*(x2-xoff))
		var ang=Math.atan2(y2-yoff, x2-xoff)
		ctx.save()
		ctx.translate(xoff, yoff)
		ctx.rotate(ang)
		ctx.drawImage(imgBond, 0, -BondWidth/2, r+5, BondWidth);
		ctx.restore()
	}
}

///
///   Clear drawing window (canvas)
///
function clear(ctx, width, height) {
	// Clear plot area
	var color="rgb(255,255,255)"
	if(drawMode==2) color="rgb(128,128,204)"
	ctx.fillStyle = color;
	ctx.fillRect(0, 0, width, height);
	//ctx.clearRect(0, 0, width, height);

	// Show border around plot area
	ctx.beginPath();
	ctx.rect(1,1,width-1,height-1);
	ctx.lineWidth = 2;
	ctx.strokeStyle = "rgb(95,95,127)";
	ctx.stroke();
	ctx.closePath();
}

///
///   Draw atoms and bonds with highlighting
///
function drawHighlight(mycanvas, sequence, AtomSize, centerx, centery) {

  // Define local variables
  var i, j, k;
  var x, y, r;
  var A;
  var maxBonds, bondto;
  var x1, y1, x2, y2;
  var molecule = Mol();
  var bonds = BondMatrix();

  // Create connection to Canvas
  var canvas = $(mycanvas);
  var ctx = canvas.getContext('2d');
  var width = canvas.width;
  var height = canvas.height;

  // Clear screen before drawing molecule
  clear(ctx, width, height);

  // Define properties of lines
  ctx.strokeStyle = "rgb(0, 0, 255)";
  ctx.lineWidth = 3.0;

  // Depth sort: Draw from back to front
  if ( typeof drawHighlight.deep == 'undefined' ) {
    drawHighlight.deep = [];
    }
  if ( drawHighlight.deep.length != molecule[0].numatoms ) {
    j=0;
    drawHighlight.deep = [];
    for (i=drawHighlight.deep.length; i<molecule[0].numatoms; i++) {
      j++;
      drawHighlight.deep[i] = [];
      drawHighlight.deep[i]["id"] = j;
      drawHighlight.deep[i]["sz"] = molecule[j].z;
      }
    }
  for (i=0; i<molecule[0].numatoms; i++) {
    j = drawHighlight.deep[i].id;
    drawHighlight.deep[i].z = molecule[j].z;
    }
  drawHighlight.deep.sort(function(a,b){
    return(b.z - a.z);
    });

  // Loop over all atoms in molecule
  for (i = 0; i < molecule[0].numatoms; i++) {
    atom = drawHighlight.deep[i].id;
    hilite = 0;
    for (ih=1; ih <= sequence[0]*3; ih++) {
      if (atom == sequence[ih])
        hilite = 1;
      }
    HighlightAtom(ctx, atom, hilite, AtomSize, centerx, centery);
    // Loop over all bonds on current atom
    for (bondto = 1; bondto <= molecule[0].numatoms; bondto++) {
      if ( bonds[atom][bondto] > 0) {
        highbond = 0;
        if (hilite > 0) {
          for (ih=1; ih < sequence[0]*3; ih++) {
            if (bondto == sequence[ih])
              highbond = 1;
            }
          }
        for (k=i+1; k < molecule[0].numatoms; k++) {
          if (drawHighlight.deep[k].id == bondto) {
            HighlightBond(ctx,atom,bondto,highbond, AtomSize, centerx, centery);
            }
          }
        }
      }
    }

  // Finished with drawHighlight routine
  }

///
///   Routine to draw a single highlighted atom
///
function HighlightAtom(ctx, atomNum, hilite, AtomSize, centerx, centery) {
  // Define local variables
  var A;
  var x, y, r;
  var xoff, yoff, roff;
  var istart, ilength, colorstr;
  var myrgb = new Array();
  var bcolor = "rgb(191,191,191)";
  var hcolor = "rgb(255,  0,  0)";
  var gcolor = "rgb(127,  0,  0)";
  var molecule = Mol();
  var activeCanvas = activeWin("");

  // Get parameters for atom
  A = molecule[atomNum].atomicnumber;
  x = molecule[0].AtomScale*molecule[atomNum].x + centerx;
  y = molecule[0].AtomScale*molecule[atomNum].y + centery;
  r = molecule[0].AtomScale*AtomSize*element(A,"radius");
  xoff = x-(0.20*r);
  yoff = y-(0.20*r);
  roff = 0.30*r;
  x = Math.floor(x+0.5);
  y = Math.floor(y+0.5);
  r = Math.floor(r+0.5);
  xoff = Math.floor(xoff+0.5);
  yoff = Math.floor(yoff+0.5);
  AtomColor = ctx.createRadialGradient(xoff, yoff, roff, x, y, r);

  alpha = 0.33;
  if (hilite > 0)
    alpha = 1.0;
  istart  = element(A,"color").indexOf("(") + 1;
  ilength = element(A,"color").indexOf(")") - istart;
  colorstr = element(A,"color").substr(istart,ilength);
  myrgb = colorstr.replace(/,/g," ").split(' ');
  hcolor = "rgba("+myrgb[0]+","+myrgb[1]+","+myrgb[2]+","+alpha+")";
  istart  = element(A,"gradient").indexOf("(") + 1;
  ilength = element(A,"gradient").indexOf(")") - istart;
  colorstr = element(A,"gradient").substr(istart,ilength);
  myrgb = colorstr.replace(/,/g," ").split(' ');
  gcolor = "rgba("+myrgb[0]+","+myrgb[1]+","+myrgb[2]+","+alpha+")";
  AtomColor.addColorStop(0, hcolor);
  AtomColor.addColorStop(1, gcolor);

  // Plot atom as a gradient shaded circle
  ctx.beginPath();
  ctx.arc(x, y, r, 0, 2*Math.PI, false);
  ctx.fillStyle = AtomColor;
  ctx.fill();
  ctx.lineWidth = 0.01;
  ctx.strokeStyle = "black";
  ctx.stroke();
  ctx.closePath();

  // Draw label if highlighted atom
  if (hilite > 0)
    atomLabel(ctx,A,atomNum,x,y);
  }

///
///   Routine to draw a single highlighted bond
///
function HighlightBond(ctx, atom1, atom2, highlite, AtomSize, centerx, centery) {
  // Define local variables
  var i, j, k;
  var x1, y1, z1, x2, y2, z2;
  var dx, dy, dz;
  var xoff, yoff;
  var dist, r;
  var AtomScale = molecule[0].AtomScale;
  var BondWidth = Math.floor(0.10*AtomScale + 0.5);
  var BondColor = "rgba(200,128,51,1.0)";
  var bcolor = "rgba(191,191,191,0.33)";
  var molecule = Mol();
  var activeCanvas = activeWin("");

  // Get parameters for atoms
  r = AtomScale * AtomSize * element(molecule[atom1].atomicnumber,"radius");
  x1 = AtomScale*molecule[atom1].x + centerx;
  y1 = AtomScale*molecule[atom1].y + centery;
  z1 = AtomScale*molecule[atom1].z;
  x2 = AtomScale*molecule[atom2].x + centerx;
  y2 = AtomScale*molecule[atom2].y + centery;
  z2 = AtomScale*molecule[atom2].z;

  // Calculate intersection of bond with sphere
  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
  xoff = x1 + r*dx/dist;
  yoff = y1 + r*dy/dist;

  // Draw bond
  ctx.beginPath();
  ctx.moveTo(xoff, yoff);
  ctx.lineTo(x2, y2);
  ctx.lineWidth = BondWidth + 2;
  ctx.strokeStyle = "rgba(0,0,0,0.33)";
  if ( highlite > 0 )
    ctx.strokeStyle = "rgba(0,0,0,1.0)";
  ctx.stroke();
  ctx.moveTo(xoff, yoff);
  ctx.lineTo(x2, y2);
  ctx.lineWidth = BondWidth;
  if ( highlite > 0 ) {
      ctx.strokeStyle = BondColor;
    } else {
      ctx.strokeStyle = bcolor;
    }
  ctx.stroke();
  ctx.closePath();
  }
// -------------------- mouse.js file --------------------

///
///   Routine to save number of atom selected with mouseDown event
///   option = "Set" or "Show"
///   value  = Number of selected atom
///        0 = MouseDown, but not on atom
///       -1 = Mouse not down
///
function mouseState(option, value) {
  // Initialize parameters as necessary
  if ( typeof mouseState.DownAtom == 'undefined' )
    mouseState.DownAtom = -1;

  // Show number of atom selected by mouseDown
  if ( option == "Show" )
    return mouseState.DownAtom;

  if ( option == "Set" )
    mouseState.DownAtom = value;

  // End of mouseState routine
  }

//#   function parameters()
//#
//#   Define constants used to control drawing.
//#   Get/set values using:
//#      var param = parameters();
//#
//#    param.mode   = Interface mode.  Set to "Draw" or "View"
//#    param.element   = Type of element to add next.  Default = "C"
//#    param.clouds = # of hybrid orbitals.  Default = 4
//#    param.bondmode   = Bond mode.  Set to "Add", "Delete", or "Rotate"
//#    param.atommode    = Add atom mode.  Set to "Add" or "Delete"
//#
function parameters() {

  // Initialize parameter array if necessary
  if ( typeof parameters.values === "undefined" )
    parameters.values = new paramObject();

  // Define default values for parameters
  if ( typeof parameters.values.mode == 'undefined' )
    parameters.values.mode = "View";
  if ( typeof parameters.values.elem == 'undefined' )
    parameters.values.elem = "C";
  if ( typeof parameters.values.clouds == 'undefined' )
    parameters.values.clouds = 4;
  if ( typeof parameters.values.add == 'undefined' )
    parameters.values.add = "Add";
  if ( typeof parameters.values.bond == 'undefined' )
    parameters.values.bond = "Add";

  // End of parameters routine
  return parameters.values;
  }

///
///   Define structure for drawing parameters
///
function paramObject() {
  this.mode = "View";
  this.element = "C";
  this.clouds = 4;
  this.bondmode = "Add";
  this.atommode = "Add";
  }

///
///   Routine to handle mouse press
///
function MouseDown(evt) {
  // Declare local variables
  var cx, cy;
  var windowName = evt.target.id;
  var activecanvas = activeWin(windowName);
  var canvas = $(activecanvas);
  var width = canvas.width;
  var height = canvas.height;

  // Map touch screen and mouse buttons
  cx = evt.pageX;
  cy = evt.pageY;
  // If touch screen
  if (evt.targetTouches) {
    if (evt.targetTouches.length == 1) {
      cx = evt.targetTouches[0].pageX;
      cy = evt.targetTouches[0].pageY;
      }
    }

  // Save selected atom
  mouseState( "Set", onAtom(cx, cy, canvas, width, height) );
  }

///
///   Routine to handle mouse up
///
function MouseUp(evt) {
  // Declare local variables
  var cx, cy;
  var mouseOnAtom, UpAtom;
  var activecanvas = activeWin("");
  var molecule = Mol();
  var canvas = $(activecanvas);
  var width = canvas.width;
  var height = canvas.height;
  var param = parameters();

  // Map touch screen and mouse buttons
  cx = evt.pageX;
  cy = evt.pageY;
  // If touch screen
  if (evt.targetTouches) {
    cx = evt.targetTouches[0].pageX;
    cy = evt.targetTouches[0].pageY;
    }

  // See which atom (if any) mouse is currently on
  mouseOnAtom = mouseState("Show");
  mouseState("Set", -1);
  Upatom = onAtom(cx, cy, canvas, width, height);

  // If 'click' moves to blank space, exit
  if ( Upatom == 0 ) {
    viewGeom(0);
    return;
    }

  // If 'click' stays on a same atom, proceed
  if ( Upatom == mouseOnAtom ) {
    // Show geometry information
    if (param.mode == "View") {
      if ( evt.shiftKey > 0) {
          molecule[Upatom].highlite = ( (molecule[Upatom].highlite+1) % 2);
          drawMolecule();

        } else {
          viewGeom(mouseOnAtom);
          }
      }
    // Add or Remove atom(s)
    if (param.mode == "Draw") {
      // Save Molecular information
      Undo("save");
      molecule[0].charge = 999;
      if (param.atommode == "Delete") {
          delAtom(Upatom);
        } else {
          newAtom(Upatom);
        }
      }
    }

  // If 'click' moves to a different atom, draw or remove bond
  if ( (Upatom != mouseOnAtom) && (param.mode == "Draw") ) {
    // Save Molecular information
    Undo("save");
    if (param.bondmode == "Add") {
      molecule[0].charge = 999;
      addBond(mouseOnAtom,Upatom);
      }
    if (param.bondmode == "Delete") {
      molecule[0].charge = 999;
      delBond(mouseOnAtom,Upatom);
      }
    if (param.bondmode == "Rotate") {
      BondRotation("Set",mouseOnAtom,Upatom);
      }
    drawMolecule();
    }

   //Finished with MouseUp routine
  }
///
///   Is mouse within radius of atom?
///   Parameters - mouse coordinates
///   Return - number of atom or zero if none
///
function onAtom(mx, my, canvas, width, height) {
  // Declare local variables
  var i, A, dx, dy;
  var radius, delta;
  var AtomSize = 0.50;
  var molecule = Mol();

  // Convert mouse position to atomic coordinate system
  mx = (mx - canvas.offsetLeft - width/2);
  my = (my - canvas.offsetTop - height/2);

  // Loop over all atoms, looking for first match
  for (i=1; i <= molecule[0].numatoms; i++) {
    A = molecule[i].atomicnumber;
    radius = molecule[0].AtomScale*AtomSize*element(A,"radius");
    radius = radius*radius;
    dx = molecule[0].AtomScale*molecule[i].x - mx;
    dy = molecule[0].AtomScale*molecule[i].y - my;
    delta = dx*dx + dy*dy;
    if (delta <= radius) {
      return i;
      }
    }
  return 0;
  }

///
///   Routine to deal with scroll wheel
///
function MouseWheel(evt) {

  // Local variable
  var molecule = Mol();

  // Stop window scrolling
  if (evt.preventDefault)
    evt.preventDefault();

  //  Firefox and Opera use detail, Chrome uses wheelDelta(?)
  var delta = evt.detail ? evt.detail*(-1) : evt.wheelDelta;
  var newScale = (delta > 0) ? 1.1 : 0.9;
  molecule[0].AtomScale = molecule[0].AtomScale*newScale;
  drawMolecule();
  return false;
  }

///
///   Alter molecule IFF mouse is pressed
///
function MouseMove(evt) {
  var newx, newy;
  var dx, dy, tx, ty, da;
  var x, y, z;
  var xc, yc;
  var cosx, cosy, cosz, sinx, siny, sinz;
  var cosA, sinA;
  var step=2;  // Smaller values make less sensitive
  var molecule = Mol();
  var param = parameters();

  // Save old (x,y) coordinates
  if ( typeof MouseMove.cx == 'undefined' )
    MouseMove.cx = 250;
  if ( typeof MouseMove.cy == 'undefined' )
    MouseMove.cy = 250;

  // Only rotate if mouse pressed, but not on an atom
  if ( mouseState("Show") == 0 ) {

    // Map touch screen and mouse buttons
    newx = evt.pageX;
    newy = evt.pageY;
    // If touch screen
    if (evt.targetTouches) {
      if (evt.targetTouches.length == 1) {
        newx = evt.targetTouches[0].pageX;
        newy = evt.targetTouches[0].pageY;
        }
      }

    // Record movement
    dx = (MouseMove.cx - newx);
    dy = (MouseMove.cy - newy);
    MouseMove.cx = newx;
    MouseMove.cy = newy;
    if ( Math.abs(dx)+Math.abs(dy) > 10)
      return;

    // Handle special case of rotation about a bond
    if ( param.bondmode == "Rotate" ) {
      da = dx + dy;
      BondRotation("Rotate",da);
      drawMolecule();
      return;
      }

    // If <Shift> key pressed, rotate around z-axis
    if ( evt.shiftKey > 0) {
      cosz = Math.cos(dx*step*Math.PI/180);
      sinz = Math.sin(dx*step*Math.PI/180);
      for (i=1; i <= molecule[0].numatoms; i++) {
        x = molecule[i].x;
        y = molecule[i].y;
        molecule[i].x =  cosz*x + sinz*y;
        molecule[i].y = -sinz*x + cosz*y;
        }
      drawMolecule();
      return;
      }

    // If <Ctrl> key pressed, translate along x- or y-axis
    if ( evt.ctrlKey > 0) {
      tx = dx*0.01;
      ty = dy*0.01;
      for (i=1; i <= molecule[0].numatoms; i++) {
        x = molecule[i].x;
        y = molecule[i].y;
        molecule[i].x = x - tx;
        molecule[i].y = y - ty;
        }
      drawMolecule();
      return;
      }

    // Rotate entire molecule
    cosx = Math.cos(dy*step*Math.PI/180);
    sinx = Math.sin(dy*step*Math.PI/180);
    cosy = Math.cos(dx*step*Math.PI/180);
    siny = Math.sin(dx*step*Math.PI/180);
    for (i=1; i <= molecule[0].numatoms; i++) {
      x = molecule[i].x;
      y = molecule[i].y;
      z = molecule[i].z;
      molecule[i].x = cosy*x + siny*z;
      molecule[i].y = -sinx*siny*x + cosx*y + sinx*cosy*z;
      molecule[i].z = -cosx*siny*x - sinx*y + cosx*cosy*z;
      }
    drawMolecule();
    }
  }

///
///   VIEW MODE: Routine to display geometry information to user
///
function viewGeom(OnAtom) {
  // Declare local variables
  var i, A1, A2, A3, A4;
  var dx, dy, dz, d;
  var dux, duy, duz, dvx, dvy, dvz, du, dv;
  var label;
  var activecanvas = activeWin("");
  var canvas = $(activecanvas);
  var ctx = canvas.getContext('2d');
  var width = canvas.width;
  var height = canvas.height;
  var molecule = Mol();

  // Initialize geom array if necessary
  if ( typeof viewGeom.geom == 'undefined' ) {
    viewGeom.geom = [0,0,0,0,0];
    }

  // If valid atom not sent, zero array
  if (OnAtom == 0) {
    for (i=0; i < 5; i++)
      viewGeom.geom[i] = 0;
    drawMolecule();
    return;
    }

  // If duplicate atom sent, zero array and return
  if (OnAtom) {
    for (i=1; i < 5; i++) {
      if (viewGeom.geom[i] == OnAtom)
        viewGeom(0);
      }
    }

  // Display appropriate information to user
  viewGeom.geom[0] = viewGeom.geom[0] + 1;
  sw = viewGeom.geom[0];
  viewGeom.geom[sw] = OnAtom;
//InfoWin("\nMouse on atom #"+OnAtom+"   Switch = "+sw);
  switch (sw*0) {
    case 1:
        A1 = molecule[viewGeom.geom[1]].atomicnumber;
        label = element(A1,"symbol") + viewGeom.geom[1].toString();
        if ( molecule[0].showcharges != 0 ) {
          label += ": Charge = ";
          if ( molecule[viewGeom.geom[1]].charge > 0 )
            label += "+";
          label += molecule[viewGeom.geom[1]].charge.toFixed(2);
          }
        drawMolecule();
        geomLabel(ctx,label,width);
        break;
    case 2:
        // Bond distance
        dx = molecule[viewGeom.geom[1]].x - molecule[viewGeom.geom[2]].x;
        dy = molecule[viewGeom.geom[1]].y - molecule[viewGeom.geom[2]].y;
        dz = molecule[viewGeom.geom[1]].z - molecule[viewGeom.geom[2]].z;
        d = Math.sqrt(dx*dx+dy*dy+dz*dz);
        d = d.toFixed(3);
        d = d.toString();
        A1 = molecule[viewGeom.geom[1]].atomicnumber;
        A2 = molecule[viewGeom.geom[2]].atomicnumber;
        label  = element(A1,"symbol") + viewGeom.geom[1].toString() + "--";
        label += element(A2,"symbol") + viewGeom.geom[2].toString();
        label += " = " + d;
        drawMolecule();
        geomLabel(ctx,label,width);
        break;
    case 3:
        // Bond angle
        A1 = molecule[viewGeom.geom[1]].atomicnumber;
        A2 = molecule[viewGeom.geom[2]].atomicnumber;
        A3 = molecule[viewGeom.geom[3]].atomicnumber;
        label  = element(A1,"symbol") + viewGeom.geom[1].toString() + "--";
        label += element(A2,"symbol") + viewGeom.geom[2].toString() + "--";
        label += element(A3,"symbol") + viewGeom.geom[3].toString();
        ang = angle(molecule,viewGeom.geom[1],viewGeom.geom[2],viewGeom.geom[3]);
        ang = ang.toFixed(1);
        ang = ang.toString();
        label += " = " + ang + "ᵒ";
        drawMolecule();
        geomLabel(ctx,label,width);
        break;
    case 4:
        // Dihedral angle
        A1 = molecule[viewGeom.geom[1]].atomicnumber;
        A2 = molecule[viewGeom.geom[2]].atomicnumber;
        A3 = molecule[viewGeom.geom[3]].atomicnumber;
        A4 = molecule[viewGeom.geom[4]].atomicnumber;
        label  = element(A1,"symbol") + viewGeom.geom[1].toString() + "--";
        label += element(A2,"symbol") + viewGeom.geom[2].toString() + "--";
        label += element(A3,"symbol") + viewGeom.geom[3].toString() + "--";
        label += element(A4,"symbol") + viewGeom.geom[4].toString();
        ang = dihedral(molecule,viewGeom.geom[1],viewGeom.geom[2],viewGeom.geom[3],viewGeom.geom[4]);
        ang = ang.toFixed(1);
        ang = ang.toString();
        label += " = " + ang;
        drawMolecule();
        geomLabel(ctx,label,width);
        for (i=0; i < 5; i++)
          viewGeom.geom[i] = 0;
        break;
      }
  }

///
///   Routine to select current element for drawing
///
function pickElem(elem) {

  var tablesize = element(1,"max");
  var i, str, mystyle;
  var mode;
  var param = parameters();
  var backcolor = "silver";
  var active = "lightskyblue";

  // If first time routine called, set defaults
  if ( typeof pickElem.metals == 'undefined' ) {
    pickElem.metals = "none";
    pickElem.myfont = "14px";
    if ($('pchooser') )
      $('pchooser').innerHTML = "Metals";
    }

  // Change display of periodic table
  if (elem == "PView") {
    if ($('pchooser') ) {
      if ( $('pchooser').innerHTML == "Metals" ) {
          pickElem.metals = "table-cell";
          pickElem.myfont = "10px";
          $('pchooser').innerHTML = "Organic";
        } else {
          pickElem.metals = "none";
          pickElem.myfont = "14px";
          $('pchooser').innerHTML = "Metals";
        }
      }
    }

  // Turn off display for metals
//  if (elem == "MOff") {
//    pickElem.metals = "none";
//    pickElem.myfont = "14px";
//    if ($('MOn') )
//      $('MOn').style.color  = "navy";
//    if ($('MOff') )
//      $('MOff').style.color = "antiquewhite";
//    if ($('pchooser') )
//      $('pchooser').innerHTML = "Metals";
//    }

  // Turn on display for metals
//  if (elem == "MOn") {
//    pickElem.metals = "table-cell";
//    pickElem.myfont = "10px";
//    if ($('pchooser') )
//      $('pchooser').innerHTML = "Organic";
//    }

  // Set display options depending on selected view
  if ($('Row1') )
      $('Row1').style.display   = pickElem.metals;
  if ($('Row2') )
      $('Row2').style.display   = pickElem.metals;
  if ($('Row3') )
      $('Row3').style.display   = pickElem.metals;
  if ($('RowLa') )
      $('RowLa').style.display  = pickElem.metals;
  if ($('RowAc') )
      $('RowAc').style.display  = pickElem.metals;
  if ($('RowLa2') )
      $('RowLa2').style.display = pickElem.metals;
  if ($('RowAc2') )
      $('RowAc2').style.display = pickElem.metals;
//  if ($('Mon') )
//      $('MOn').style.fontSize   = pickElem.myfont;
//  if ($('Moff') )
//      $('MOff').style.fontSize  = pickElem.myfont;

  // Loop over all elements in table
  for (i=1; i<tablesize; i++) {
    delete mystyle;
    // Process Main Group elements
    str = "m" + element(i,"symbol");
    if ($(str)) {
      mystyle = $(str).style;
      }
    // Process metallic elements (not Main Group)
    str = "p" + element(i,"symbol");
    if ($(str)) {
      mystyle = $(str).style;
      mystyle.display = pickElem.metals;
      }
    // Set display parameters for all elements
    if ( mystyle ) {
      mystyle.backgroundColor = backcolor;
      mystyle.fontSize = pickElem.myfont;
      if (element(i,"symbol") == elem)
        mystyle.backgroundColor = active;
      }
    }

  // Set display parameters for all buttons
  buttonColor("delbtn",0);
  buttonColor("delbond",0);
  buttonColor("rotbond",0);
  if (elem == "Delete") {
    buttonColor("delbtn",1);
    param.atommode="Delete";
    return;
    }
  if (elem == "DelBond") {
    buttonColor("delbond",1);
    param.bondmode="Delete";
    return;
    }
  if (elem == "RotateBond") {
    if (  param.bondmode != "Rotate" ) {
        buttonColor("rotbond",1);
        param.bondmode="Rotate";
      } else {
        param.bondmode="";
        BondRotation("Clear");
      }
    return;
    }

  // Not a button, so set element
  param.element = elem;
  param.atommode = "Add";
  param.bondmode = "Add";

  // Finished with pickElem routine
  }

//#   function viewmode()
//#
//#   Routine to enable display of "view mode" interface.
//#
function viewmode() {
  var param = parameters();
  if ( $("drawDiv") )
    $("drawDiv").style.display = "none";
  if ( $("viewDiv") )
    $("viewDiv").style.display = "inline";
  buttonColor("modeDraw",0);
  buttonColor("modeView",1);
  param.mode = "View";
  pickElem("ClearBond");

  // Show molecular formula
  formula();
  //showCoord(1);
  }

///
///   Routine to allow rotation about a bond
///
///   mode: "Set" sets atoms, "Rotate" performs rotation, or "Clear"
///   values: Two atoms of bond  OR  single rotation angle
///
function BondRotation(mode,value1,value2) {

  // Define local variables
  var i, pos;
  var x, y, xc, yc;
  var cosA, sinA;
  var step=1.0;  // Larger values make more sensitive
  var molecule = Mol();

  // Save two atoms used to define bond rotation
  if ( typeof BondRotation.RotateAtom1 == 'undefined' )
    BondRotation.RotateAtom1 = 0;
  if ( typeof BondRotation.RotateAtom2 == 'undefined' )
    BondRotation.RotateAtom2 = 0;

  // Save list of atoms to rotate if rotation about a bond
  if ( typeof BondRotation.RotateList == 'undefined' ) {
    BondRotation.RotateList = new Array();
    BondRotation.RotateList[0] = 0;
    }

  // If mode is "Clear", reset variables
  if ( mode == "Clear" ) {
    BondRotation.RotateAtom1 = 0;
    BondRotation.RotateAtom2 = 0;
    BondRotation.RotateList[0] = 0;
    return;
    }

  // If mode is "Set", align molecule along bond
  if ( mode == "Set" ) {
    BondRotation.RotateAtom1 = value1;
    BondRotation.RotateAtom2 = value2;
    bondAlign(value1,value2);
    return;
    }

  // If mode is "Rotate", perform rotation
  if ( mode == "Rotate" ) {
    atom1 = BondRotation.RotateAtom1;
    atom2 = BondRotation.RotateAtom2;
    cosA = Math.cos(value1*step*Math.PI/360);
    sinA = Math.sin(value1*step*Math.PI/360);
    if (BondRotation.RotateList[0] == 0) {
      BondRotation.RotateList[0] = 1;
      BondRotation.RotateList[1] = BondRotation.RotateAtom2;
      BondRotation.RotateList = rotlist(BondRotation.RotateList, atom1, atom2);
      }

InfoWin("\n--- Bond rotation ---\n");
for(i=1;i<=BondRotation.RotateList[0];i++){
 j=BondRotation.RotateList[i];
 InfoWin(element(molecule[j].atomicnumber,"symbol")+j+"\n");
 }

    xc = molecule[atom1].x;
    yc = molecule[atom1].y;
    for (pos=1; pos <= BondRotation.RotateList[0]; pos++) {
      i = BondRotation.RotateList[pos];
      x = molecule[i].x - xc;
      y = molecule[i].y - yc;
      molecule[i].x =  cosA*x + sinA*y + xc;
      molecule[i].y = -sinA*x + cosA*y + yc;
      }
    }

  // End of BondRotation routine
  }

///
///   This routine aligns molecule so that bond is parallel to z-axis
///
function bondAlign(atom1, atom2) {
  var i;
  var x, y, z;
  var dx, dy, dz;
  var ctx, cty, stx, sty;
  var molecule = Mol();

  // Determine elements of rotation matrix used to align bond
  dx = molecule[atom2].x - molecule[atom1].x;
  dy = molecule[atom2].y - molecule[atom1].y;
  dz = molecule[atom2].z - molecule[atom1].z;
  if ( (dy == 0) && (dz == 0) ) {
      ctx = 1.0;
      stx = 0.0;
    } else {
      ctx = Math.sqrt( dz*dz / ((dy*dy)+(dz*dz)) );
      stx = Math.sqrt( 1 - (ctx*ctx) );
    }
  vsign = dy*dz;
  if (vsign < 0.0)
    stx = -stx;
  z = stx*dy + ctx*dz;
  if ( (dx == 0) && (z == 0) ) {
      cty = 1.0;
      sty = 0.0;
    } else {
      cty = Math.sqrt( z*z / ((dx*dx)+(z*z)) );
      sty = Math.sqrt( 1 - (cty*cty) );
    }
  vsign = dx*z;
  if (vsign < 0.0)
    sty = -sty;

  // Rotate molecule
  for (i=1; i <= molecule[0].numatoms; i++) {
    x = molecule[i].x;
    y = molecule[i].y;
    z = molecule[i].z;
    molecule[i].x = (x*cty) - (y*stx*sty) - (z*ctx*sty);
    molecule[i].y =           (y*ctx)     - (z*stx);
    molecule[i].z = (x*sty) + (y*stx*cty) + (z*ctx*cty);
    }

  // Make sure second atom in front
  if ( molecule[atom1].z < molecule[atom2].z ) {
    for (i=1; i <= molecule[0].numatoms; i++) {
      molecule[i].x = -molecule[i].x;
      molecule[i].z = -molecule[i].z;
      }
    }

  // End of rotateBond routine
  }

///
///   Routine to determine which atoms bonded to atom2,
///   after but not including atom1
///
function rotlist(RList, atom1, atom2) {
  // Declare local variables
  var i, j;
  var atm, found;
  var molecule = Mol();
  var bonds = BondMatrix();

  // Loop over all bonds on atom2
  for (atm=1; atm <= molecule[0].numatoms; atm++) {
    if ( bonds[atom2][atm] > 0) {
      found = 0;
      if ( atm == atom1 )
        found = 1;
      for (j=1; j <= RList[0]; j++) {
        if ( atm == RList[j] ) {
          found = 1;
          j = RList[0] + 1;
          }
        }
      if (found == 0) {
        RList[0]++;
        RList[RList[0]] = atm;
        rotlist(RList, atom2, atm);
        }
      }
    }

  // End of rotlist routine
  return RList;
  }

// -------------------- End of mouse.js file --------------------
//
// -------------------- Section to handle Reading Files --------------------

///
///   Routine to get filename of local input file from user
///
function MyFileReader(evt) {
  var validExtension;
  var localfile = evt.target.files[0];
  // Call routine to read and process data
  if (localfile) {
      var fr = new FileReader();
      fr.onload = function(fh) {
        validExtension = 0;
        var ext = extension(localfile.name);
        // Read contents and split data into lines
        var contents = fh.target.result;
        var lines = contents.split("\n");
        if ( ext == "pdb" ) {
          validExtension = 1;
          readPDBfile(lines);
          }
        if ( ext == "xyz" ) {
          validExtension = 1;
          readXYZfile(lines);
          }
        if ( ext == "mol" ) {
          validExtension = 1;
          readMOLfile(localfile,lines);
          }
        if ( ext == "inp" ) {
          validExtension = 1;
          readINPfile(localfile,lines);
          }
        if ( validExtension == 0 ) {
          InfoWin("*** ERROR: Invalid file type for "+localfile.name+". ***",1);
          return;
          }
        centerMolecule();
        drawMolecule();
        }
      fr.readAsText(localfile);
    } else {
      alert("Failed to load file");
    }
  }

///
///   Routine to determine extension for filename
///
function extension(myfilename) {
  var extension = myfilename.substring(myfilename.lastIndexOf('.')+1);
  return extension.toLowerCase();
  }

//#
//#   function InfoWin(mytext,mode)
//#
//#   Routine to write text to information window (textarea with an id of "information").
//#
//#   Parameters:
//#     mytext:  String to write to output window
//#     mode:    If >0, then clear text window
//#
function InfoWin(mytext,mode) {
  var outstring = mytext || "";
  var clearwin = mode || 0;

  if ( ! $("information") )
    return;
  if ( clearwin > 0 )
    $("information").value  = outstring;
  else
    $("information").value += outstring;
  }

//#   function loadMolecule()
//#
//#   Routine to load molecular information for methane molecule.
//#
function loadMolecule() {

  // Add atomic coordinates to molecule array
  addAtom(6,  0.000,  0.000,  0.000);
  addAtom(1,  0.874,  0.618,  0.000);
  addAtom(1, -0.874,  0.618,  0.000);
  addAtom(1,  0.000, -0.618,  0.874);
  addAtom(1,  0.000, -0.618, -0.874);
  // Add bonds to bond array
  addBond(1, 2);
  addBond(1, 3);
  addBond(1, 4);
  addBond(1, 5);
  }

//#   function resetView()
//#
//#   Routine to reset view. Center, rescale, and remove highlights.
//#
function resetView() {

  // Define local variables
  var i;
  var molecule = Mol();

  // Remove highlights and show all atoms
  for ( i=1; i <= molecule[0].numatoms; i++ ) {
    molecule[i].highlite = 0;
    molecule[i].hide = 0;
    }

  // Clear Information Window
  InfoWin("",1);

  // Formula, center, scale, and draw
  formula();
  centerMolecule();
  drawMolecule();
  }

//#   function buttonColor(button, mode)
//#
//#   Set color of buttons
//#
//#   Parameters
//#     button - string containing ID name for button
//#     mode   - 0=inactive, 1=active
//#
function buttonColor(button, mode) {

  mode = mode || 0;
  var buttoncolor = "silver";
  if ( mode )
    buttoncolor = "lightskyblue";
  if ( $(button) )
    $(button).style.backgroundColor = buttoncolor;
  }
