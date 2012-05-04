
/*!
 * simplexn
 * dimension-independent geometric kernel based on simplicial complex
 * Copyright (c) 2011 cvd-lab <cvd-lab@email.com> (https://github.com/cvd-lab/)
 * MIT License
 */

!(function (exports) {

  /**
   * Variables.
   */

  var cos = Math.cos;
  var sin = Math.sin;
  var round = Math.round;
  var min = Math.min;
  var abs = Math.abs;
  var pi = Math.PI;
  var random = Math.random;
  var floor = Math.floor;

  /**
   * Library namespace.
   */

  var simplexn = exports.simplexn = {};

  /**
   * Library version.
   */

  simplexn.version = '0.1.5';

  /**
   * utils namespace
   * @api private
   */

  simplexn._utils = {};

  /**
   * _flat
   * Return a flat version of the given array of arrays.
   * 
   * @param {Array} arrays
   * @return {Array} array
   * @api private
   */

  var _flat =
  simplexn._utils._flat = function (arrays) {
    var res = [];

    arrays.forEach(function (item) {
      res = res.concat(item);
    });

    return res;
  };

  /**
   * _repeat
   * Return an array made by n times value item.
   * 
   * @param {Number|Boolean|String} value
   * @param {Number} n
   * @return {Array} array
   * @api private
   */

  var _repeat = 
  simplexn._utils._repeat = function (value, n) {
    var res = [];

    while (n--) res.push(value);

    return res;
  };

  /**
   * _swap
   * Swap i1 to i2 indexed items in array.
   * 
   * @param {Array|BufferArray} array
   * @param {Number} i1
   * @param {Number} i2
   * @api private
   */

  var _swap = 
  simplexn._utils._swap = function (array, i1, i2) {
    var tmp = array[i1];
    array[i1] = array[i2];
    array[i2] = tmp;
  };

  /**
   * _quickSort
   * Quick sort algorithm.
   *
   * @param {Array|BufferArray} array to sort
   * @api private
   */

  var _partition = 
  simplexn._utils._partition = function (array, begin, end, pivot) {
    var piv = array[pivot];
    var store = begin;
    var ix;

    _swap(array, pivot, end - 1);
    for (ix = begin; ix < end - 1; ++ix) {
      if (array[ix] <= piv) {
        _swap(array, store, ix);
        ++store;
      }
    }
    _swap(array, end - 1, store);

    return store;
  };

  var _qsort = 
  simplexn._utils._qsort = function (array, begin, end) {
    if (end - 1 > begin) {
      var pivot = begin + floor(random() * (end - begin));

      pivot = _partition(array, begin, end, pivot);

      _qsort(array, begin, pivot);
      _qsort(array, pivot + 1, end);
    }
  };

  var _quickSort = 
  simplexn._utils._quickSort = function (array) {
    _qsort(array, 0, array.length);
  };

  /**
   * _areEqual
   * 
   * @param {Array|Float32Array|Uint32Array} a1
   * @param {Array|Float32Array|Uint32Array} a2
   * @return {Boolean} true if each item of a1 is === to correspond element of a2
   * @api private
   */

  var _areEqual = 
  simplexn._utils._areEqual = function (a1, a2) {
    var a1Len = a1.length;
    var a2Len = a2.length;
    var i;

    if (a1Len !== a2Len) {
      return false;
    }

    for (i = 0; i < a1Len; i++) {
      if (a1[i] !== a2[i]) {
        return false;
      }
    }
 
    return true;
  };

  /**
   * vector operations namespace
   * @api public
   */

  simplexn.vector = {};

  /**
   * add
   * 
   * @param {Array|Float32Array|Uint32Array} v1
   * @param {Array|Float32Array|Uint32Array} v2
   * @return {Array|Float32Array|Uint32Array}
   * @api public
   */

  var vectorAdd =
  simplexn.vector.add = function (v1, v2) {
    var rn = v1.length;
    var res = new v1.constructor(rn);
    var i;

    for (var i = 0; i < rn; i += 1) {
      res[i] = v1[i] + v2[i];
    };

    return res;
  };

  /**
   * sub
   * 
   * @param {Array|Float32Array|Uint32Array} v1
   * @param {Array|Float32Array|Uint32Array} v2
   * @return {Array|Float32Array|Uint32Array}
   * @api public
   */

  var vectorSub =
  simplexn.vector.sub = function (v1, v2) {
    var rn = v1.length;
    var res = new v1.constructor(rn);
    var i;

    for (var i = 0; i < rn; i += 1) {
      res[i] = v1[i] - v2[i];
    };

    return res;
  };

  /**
   * mul
   * 
   * @param {Array|Float32Array|Uint32Array} v1
   * @param {Array|Float32Array|Uint32Array} v2
   * @return {Array|Float32Array|Uint32Array}
   * @api public
   */

  var vectorMul =
  simplexn.vector.mul = function (v1, v2) {
    var rn = v1.length;
    var res = new v1.constructor(rn);
    var i;

    for (var i = 0; i < rn; i += 1) {
      res[i] = v1[i] * (v2[i] || 1);
    };

    return res;
  };

  /**
   * scalarMul
   * 
   * @param {Number} scalar
   * @param {Array|Float32Array|Uint32Array} v
   * @return {Array|Float32Array|Uint32Array}
   * @api public
   */

  var vectorScalarMul =
  simplexn.vector.scalarMul = function (scalar, v) {
    var rn = v.length;
    var res = new v.constructor(rn);
    var i;

    for (var i = 0; i < rn; i += 1) {
      res[i] = scalar * v[i];
    };

    return res;
  };

  /**
   * scalarDiv
   * 
   * @param {Number} scalar
   * @param {Array|Float32Array|Uint32Array} v
   * @return {Array|Float32Array|Uint32Array}
   * @api public
   */

  var vectorScalarDiv =
  simplexn.vector.scalarDiv = function (scalar, v) {
    var rn = v.length;
    var res = new v.constructor(rn);
    var i;

    for (var i = 0; i < rn; i += 1) {
      res[i] = v[i] / scalar;
    };

    return res;
  };

  /**
   * average
   * 
   * @param {Array} vectors
   * @return {Array|Float32Array|Uint32Array}
   * @api public
   */

  var vectorAvg =
  simplexn.vector.avg = function (vectors) {
    var vectors = vectors || [[]];
    var length = vectors.length;
    var rn = vectors[0].length;
    var res = new vectors[0].constructor(rn);
    var i;

    for (i = 0; i < rn; i += 1) {
      res[i] = 0;
    }

    res = vectors.reduce(vectorAdd, res);
    res = vectorScalarDiv(length, res);

    return res;
  };

  /**
   * matrix operations namespace
   * @api public
   */

  simplexn.matrix = {};

  /**
   * identity
   * 
   * @param {Number} dim
   * @api public
   */

  var matrixIdentity =
  simplexn.matrix.identity = function (dim) {
    var matrix = new Array(dim);
    var i, j;

    for (i = 0; i < dim; i += 1) {
      matrix[i] = new Array(dim);
      for(j = 0; j < dim; j += 1) {
        matrix[i][j] = (j === i) ? 1 : 0;
      }
    }

    return matrix;
  };

  /**
   * PointSet
   * 
   * @constructor
   * @param {Array|Number} points or number of point to initialize;
   * @rn {Number} [rn=points[0].length] points dimension;
   * @api public
   */

  var PointSet = 
  simplexn.PointSet = function (points, rn) {
    points = points || [[]];
    if (typeof points === 'number') {
      this.rn = rn;
      points = points * rn;
    } else {
      this.rn  = points[0].length;
      points = _flat(points);
    }
    this.points = new Float32Array(points);
  };

  /**
   * size
   * 
   * @property
   * @api public
   */

  simplexn.PointSet.prototype.__defineGetter__('size', function () {
    return this.points.length / this.rn;
  });

  /** 
   * clone
   * 
   * @return {PointSet} clone
   * @api public
   */

  simplexn.PointSet.prototype.clone = function () {
    var clone = new PointSet();
    clone.size = this.size;
    clone.rn = this.rn;
    clone.points = new Float32Array(this.points);
    return clone;
  };

  /**
   * equals 
   * 
   * @param {simplexn.PointSet} pointSet
   * @return {Boolean} true if this is equals to the given point set, false otherwise.
   * @api public
   */

  simplexn.PointSet.prototype.equals = function (other) {
    if (this.rn !== other.rn || this.size !== other.size) return false;
    for (var i = 0, l = this.points.length; i < l; i += 1) {
      if (this.points[i] !== other.points[i]) {
        return false;
      }
    }
    return true;
  };


  /**
   * get
   * 
   * @param {Number} index
   * @return {Float32Array} the indexed point
   * @api public
   */

  simplexn.PointSet.prototype.get = function (index) {
    var rn = this.rn;
    var begin = index * rn;
    var end = begin + rn;

    return this.points.subarray(begin, end);
  };

  /**
   * set
   * 
   * @param {Array|Float32Array} points
   * @param {Number} [index=0]
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.set = function (point, index) {
    point = point || 0;
    this.points.set(point, index * this.rn);
    return this;
  };

  /**
   * forEach
   * 
   * @param {Function} iterator
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.forEach = function (iterator) {
    var points = this.points;
    var length = points.length;
    var rn = this.rn;
    var i, j;

    for (i = j = 0; i < length; i += rn, j += 1) {
      iterator(points.subarray(i, i + rn), j);
    }

    return this;
  };

  /**
   * map
   * 
   * @param {Function} mapping
   * @return {simplexn.PointSet} a new point set
   * @api public
   */

  simplexn.PointSet.prototype.map = function (mapping) {
    var points = this.points;
    var oldRn = this.rn;
    var size = this.size;
    var mappedPoints0 = mapping(points.subarray(0,oldRn));
    var newRn = mappedPoints0.length;
    var newPoints = new Float32Array(size * newRn);
    var i, j;

    newPoints.set(mappedPoints0);

    for (i = oldRn, j = 1; j < size; i += oldRn, j += 1) {
      newPoints.set(mapping(points.subarray(i, i + oldRn), j), j * newRn);
    }

    this.points = newPoints;
    this.rn = newRn;
    return this;
  };

  /**
   * filter
   * 
   * @param {Function} iterator
   *   
   * @return {Float32Array} new filtered PointSet
   * @api public
   */

  simplexn.PointSet.prototype.filter = function (iterator) {
    var points = this.points;
    var length = points.length;
    var filtered = new Float32Array(length);
    var rn = this.rn;
    var i, j, k;
    var point;
    var pointset;

    for (i = j = k = 0; i < length; i += rn, j += 1) {
      point = points.subarray(i, i + rn);
      if (iterator(point, j)) {
        filtered.set(point, k);
        k += rn;
      }
    }

    filtered = filtered.subarray(0, k);
    pointset = new PointSet();
    pointset.points = filtered;
    pointset.rn = rn;
    pointset.size = k / rn;
    
    return pointset;
  };

  /**
   * merge
   * Filter duplicated and overlapped vertices 
   * according to precision parameter (10e-4 by default).
   * 
   * @param {Number} [precision = 10e-4] 
   * @return {Float32Array} inidices mapping changes
   * @api public
   */

  simplexn.PointSet.prototype.merge = function (precision) {
    var precision = precision || 1e-4;
    var points = this.points;
    var length = points.length;
    var rn = this.rn;
    var size = this.size;
    var indices = new Uint32Array(size);
    var merged = new Float32Array(length);
    var usedIndices = 0;
    var usedCoords = 0;
    var vertexAdded;
    var equals;
    var i, j, k;

    for (i = 0; i < length; i += rn) {
      vertexAdded = false;
      for (j = 0; j < usedCoords && !vertexAdded; j += rn) {
        equals = true;
        for (k = 0; k < rn; k += 1) {
          points[i+k] = round(points[i+k] / precision) * precision;
          equals &= points[i+k] === merged[j+k];
        }
        vertexAdded |= equals; 
      }
      indices[i/rn] = !vertexAdded ? usedIndices : j/rn-1;
      if (!vertexAdded) {
        for (k = 0; k < rn; k += 1) {
          merged[usedCoords+k] = points[i+k];
        }
        usedIndices += 1;
        usedCoords = usedIndices*rn;
      }
    }

    this.points = merged.subarray(0, usedCoords);

    return indices;
  };

  /**
   * rotate
   * a 3d rotation
   * 
   * @param {Array|Uint32Array} dims
   * @param {Number} angle
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.rotate = function (dims, angle) {
    var dims = dims[0] > dims[1] ? [dims[1], dims[0]] : dims;
    var points = this.points;
    var length = points.length;
    var rn = this.rn;
    var cos_a = cos(angle);
    var sin_a = sin(angle);
    var r_ii = cos_a;
    var r_ij = -sin_a;
    var r_ji = sin_a;
    var r_jj = cos_a;
    var d_i = dims[0];
    var d_j = dims[1];
    var v_i;
    var v_j;
    var i, j, k;

    if ((dims[0] + dims[1]) % 2 == 0) {
      r_ij *= -1;
      r_ji *= -1;
    }

    for (k = 0, i = d_i, j = d_j; k < length; k += rn, i = k + d_i, j = k + d_j) {
      v_i = points[i];
      v_j = points[j];
      points[i] = v_i * r_ii + v_j * r_ij;
      points[j] = v_i * r_ji + v_j * r_jj;
    }
    
    return this;
  };

  /**
   * scale
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.scale = function (dims, values) {
    var points = this.points 
    var length = points.length;
    var dimsLength = dims.length;
    var rn = this.rn;
    var i, j;

    for (i = 0; i < length; i += rn) {
      for (j = 0; j < dimsLength; j += 1) {
        points[i+dims[j]] *= values[j]; 
      }
    }

    return this;
  };

  /**
   * translate
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.translate = function (dims, values) {
    var rn = this.rn;
    var maxDim = Math.max.apply(null, dims.concat(rn - 1));
    this.embed(maxDim + 1);

    var points = this.points 
    var length = points.length;
    var dimsLength = dims.length;
    var i, j;

    for (i = 0; i < length; i += rn) {
      for (j = 0; j < dimsLength; j += 1) {
        points[i+dims[j]] += values[j]; 
      }
    }
    return this;
  };

  /**
   * transform
   * 
   * @param {Array|Float32Array} matrix
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.transform = function (matrix) {
    // body...
    
    return this;
  };

  /**
   * embed
   * 
   * @param {Number} dim
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.embed = function (dim) {
    var dim = dim || this.rn + 1
    var rn = this.rn;
    var minDim = Math.min(rn, dim);
    var oldPoints = this.points;
    var oldLength = oldPoints.length;
    var length = oldLength / rn * dim;
    var points = new Float32Array(length);
    var i, j, k;

    for (i = 0, j = 0; i < oldLength; i += rn, j += dim) {
      for (k = 0; k < minDim; k += 1) {
        points[j + k] = oldPoints[i + k];
      }
    }

    this.points = points;
    this.rn = dim;
    
    return this;
  };

  /**
   * prod
   * Execute product of this pointset pointset.
   *
   * @param {simplexn.PointSet} pointset
   * @return {simplexn.PointSet} this for chaining
   * @api public
   */

  simplexn.PointSet.prototype.prod = function(pointset) {
    var size = this.size;
    var rn = this.rn
    var pointsetSize = pointset.size;
    var pointsetRn = pointset.rn;
    var newSize = size * pointsetSize;
    var newRn = rn + pointsetRn;
    var newLength = newSize * newRn;
    var newPoints = new Float32Array(newLength);
    var newPoint, point1, point2;
    var i, j;
    var n = 0;

    for (j = 0; j < pointsetSize; j += 1) {
      point2 = pointset.get(j);
      for (i = 0; i < size; i += 1) {
        newPoint = new Float32Array(newRn);
        point1 = this.get(i);
        newPoint.set(point1);
        newPoint.set(point2, rn);
        newPoints.set(newPoint, newRn*n++);
      }
    }

    this.points = newPoints;
    this.rn = newRn;

    return this;
  };

  /**
   * Topology
   * 
   * @constructor
   * @param {Array|Uint32Array} complex
   * @param {Number} [dim=complex[0].length - 1]
   * @api public
   */

  var Topology = 
  simplexn.Topology = function (complex, dim) {
    this._computeTopology(complex, dim);
  }

  /**
   * _computeTopologyogy
   * 
   * @param {Array|Uint32Array} complex
   * @param {Number} [dim=complex[0].length - 1]
   * @api private
   */

  simplexn.Topology.prototype._computeTopology = function (complex, dim) {
    var complex = complex || [[]];
    var dim;
    var complexes = new Array();
    var complexTemp, complexNext;
    var complexNextLength;
    var complexLength;
    var cellDim;
    var d, c, i, j, k;
    var exchange1, exchange2;

    complex = complex.length > 0 ? complex : [[]];
    dim = dim || complex[0].length - 1;
    complex = complex instanceof Array ? _flat(complex) : complex;
    if (dim >= 0) { complexes[0] = new Uint32Array(); }
    if (dim === 0 || dim >= 1) { complexes[dim] = new Uint32Array(complex); }

    for (d = dim; d > 1; d -= 1) {
      complexNext = complexes[d];
      complexNextLength = complexNext.length;
      cellDim = d + 1;
      complexLength = complexNextLength / cellDim;
      complexTemp = new Uint32Array(cellDim * complexLength * d);
      complexes[d-1] = complexTemp;
      k = 0;
      for (c = 0; c < complexNextLength; c += cellDim) {
        for (i = 0; i < cellDim; i++) {
          for (j = 0; j < cellDim; j++) {
            if (i != j) {
              complexTemp[k] = complexNext[c+j];
              k++;
            }
          }
          if (i & 1) { // is odd
            // exchange1 = k - cellDim + 1;
            exchange1 = k - 2;
            exchange2 = k - 1;
            complexTemp[exchange1] ^= complexTemp[exchange2];
            complexTemp[exchange2] ^= complexTemp[exchange1];
            complexTemp[exchange1] ^= complexTemp[exchange2];
          }
        }
      }
    }

    this.complexes = complexes;
  };

  /**
   * dim
   * 
   * @property
   * @api public
   */

  simplexn.Topology.prototype.__defineGetter__('dim', function () {
    return this.complexes.length - 1;
  });

  /**
   * maxCells
   * 
   * @property
   * @api public
   */

  simplexn.Topology.prototype.__defineGetter__('maxCells', function () {
    return this.complexes[this.dim];
  });

  /**
   * equals 
   * 
   * @param {simplexn.Topology} topology
   * @return {Boolean} true if this is equals to the given topology, false otherwise.
   * @api public
   */

  simplexn.Topology.prototype.equals = function (other) {
    var complexes1 = this.complexes;
    var complexes2 = other.complexes;
    var dim1 = this.dim;
    var dim2 = other.dim;
    var i;

    if (dim1 !== dim2) return false;
    for (i = 0; i < dim1; i += 1) {
      if (!_areEqual(complexes1[i], complexes2[i])) return false;
    }
    return true;
  };

  /**
   * remap
   *
   * Remap topology by given mapping array.
   *
   * @param {Array} mapping
   * @return {simplexn.Topology} cloned topology
   * @api public
   */

  simplexn.Topology.prototype.remap = function (mapping) {
    var length;
    var i;

    this.complexes.forEach(function (complex) {
      length = complex.length;
      for (var i = 0; i < length; i += 1) {
        complex[i] = mapping[complex[i]];
      }
    });

    return this;
  };

  /**
   * invert
   * Invert orientation of all of the cells in topology
   * 
   * @return {simplexn.Topology} this for chaining 
   * @api public
   */

  simplexn.Topology.prototype.invert = function () {
    var dim = this.dim;
    var complex = this.complexes[dim];
    var length = complex.length;
    var cellSize = dim + 1;
    var cells = [];
    var cell;
    var swap;
    var i, j;

    for (i = 0; i < length; i += cellSize) {
      cell = [];
      for (j = 0; j < cellSize; j += 1) {
        cell.push(complex[i+j]);
      }
      swap = cell[0];
      cell[0] = cell[dim];
      cell[dim] = swap; 
      cells.push(cell);
    }

    this._computeTopology(cells);
  };

  /**
   * clone
   * 
   * @return {simplexn.Topology} cloned topology
   * @api public
   */

  simplexn.Topology.prototype.clone = function () {
    var clone = new Topology();
    var dim = this.dim;
    var complexes = new Array();
    var i;

    this.complexes.forEach(function (complex, i) {
      complexes[i] = new Uint32Array(complex);
    });

    clone.complexes = complexes;

    return clone;
  };

  /**
   * cells0d
   * 
   * @return {Uint32Array} 0-dimension cells
   * @api private
   */
   
  simplexn.Topology.prototype.cells0d = function () {
    if (this.dim === 0) { return this.complexes[0]; }

    var complexes = this.complexes || [[]];
    var cells1d = complexes[1] || [];
    var length = cells1d.length;
    var cells0d = new Uint32Array(length);
    var i, j;
    var k = 0;
    var found;

    for (i = 0; i < length; i += 1) {
      found = false;
      for (j = 0; j < k && !found; j += 1) {
        found |= cells1d[i] === cells0d[j];
      }
      if (!found) {
        cells0d[k++] = [cells1d[i]];
      }
    }

    return cells0d.subarray(0,k);
  }

  /**
   * skeleton
   * 
   * @param {Number} ord skeleton order
   * @return {simplexn.Topology} topology of the skeleton
   * @api public
   */

  simplexn.Topology.prototype.skeleton = function (ord) {
    var dim = this.dim
    var ord = ord === undefined ? dim - 1 : ord;
    var out = dim - ord;

    if (ord === 0) {
      this.complexes = [this.cells0d()];
    } else {
      while (out--) this.complexes.pop();
    }

    return this;
  };


  /**
   * boundary
   * 
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.Topology.prototype.boundary = function () {
    var dim = this.dim - 1;

    this.skeleton(dim);
    
    var complexes = this.complexes;
    var cellLength = dim + 1;
    var cells = complexes[dim];
    var cellsLength = cells.length;
    var cellsSize = cellsLength / cellLength;
    var sortedCells = new Uint32Array(cells);
    var notBoundaryCells = new Uint8Array(cellsSize);
    var boundary;
    var boundarySize = cellsSize;
    var cell;
    var equal;
    var i, j, b, c;

    for (i = 0; i < cellsLength; i += cellLength) {
      _quickSort(sortedCells.subarray(i, i + cellLength));
    }
    for (c = 0; c < cellsSize; c += 1) {
      cell = sortedCells.subarray(c * cellLength, c * cellLength + cellLength);
      if (!notBoundaryCells[c]) {
        for (i = c + 1; i < cellsSize; i += 1) {
          equal = true;
          for (j = 0; j < cellLength && equal; j += 1) {
            equal &= sortedCells[i*cellLength+j] === cell[j]; 
          }
          notBoundaryCells[c] |= equal;
          notBoundaryCells[i] |= equal;
        }
      }
    }
    for (c = 0; c < cellsSize; c += 1) {
      boundarySize -= notBoundaryCells[c];
    }
    boundary = new Uint32Array(boundarySize * cellLength);
    for (c = 0, b = 0; c < cellsSize; c += 1) {
      if (!notBoundaryCells[c]) {
        for (i = 0; i < cellLength; i += 1) {
          boundary[b++] = cells[c*cellLength+i];
        }
      }
    }

    this._computeTopology(boundary, dim);
    return this;
  };

  /**
   * SimplicialComplex
   * 
   * @constructor
   * @param {Array|Float32Array} points
   * @faces {Array|Uint32Array} complex
   * @api public
   */

  var SimplicialComplex =
  simplexn.SimplicialComplex = function (points, complex) {
    var points = points || [[]];
    var complex = complex || [[]];

    this.pointset = new PointSet(points);
    this.topology = new Topology(complex);
  };

  /**
   * rn
   * 
   * @property
   * @api public
   */

  simplexn.SimplicialComplex.prototype.__defineGetter__('rn', function () {
    return this.pointset.rn;
  });

  /**
   * size
   * 
   * @property
   * @api public
   */

  simplexn.SimplicialComplex.prototype.__defineGetter__('size', function () {
    return this.pointset.size;
  });

  /**
   * dim
   * 
   * @property
   * @api public
   */

  simplexn.SimplicialComplex.prototype.__defineGetter__('dim', function () {
    return this.topology.dim;
  });

  /**
   * rotate
   * 
   * @param {Array|Uint32Array} dims
   * @param {Number} angle
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.rotate = function (dims, angle) {
    this.pointset.rotate(dims, angle);
    return this;
  };

  /**
   * scale
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.scale = function (dims, values) {
    var invert;

    this.pointset.scale(dims, values);

    invert = values.reduce(function (v1, v2) {
      return v1 * v2;
    });

    if (invert < 0) {
      this.topology.invert();
    };
    
    return this;
  };

  /**
   * translate
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.translate = function (dims, values) {
    this.pointset.translate(dims, values);
    return this;
  };

  /**
   * transform
   * 
   * @param {Array|Float32Array} matrix
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.transform = function (matrix) {
    this.pointset.transform(matrix);
    return this;
  };

  /**
   * embed
   * 
   * @param {Number} dim
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.embed = function (dim) {
    this.pointset.embed(dim);
    return this;
  };

  /**
   * merge
   * 
   * @precision {Number} [precision=1e-4]
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.merge = function (precision) {
    var precision = precision || 1e-4;
    var mapping = this.pointset.merge(precision);
    this.topology.remap(mapping);
    return this;
  };

  /**
   * clone
   * 
   * @return {simplexn.SimplicialComplex} cloned SimplicialComplex
   * @api public
   */

  simplexn.SimplicialComplex.prototype.clone = function () {
    var clone = new SimplicialComplex();
    clone.pointset = this.pointset.clone();
    clone.topology = this.topology.clone();
    return clone;
  };

  /**
   * map
   * 
   * @param {Function} mapping
   * @param {Boolean|Number} merge
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.map = function (mapping, merge) {
    var precision = typeof merge === 'boolean' ? undefined : merge;
    this.pointset.map(mapping);
    if (merge) this.merge(precision);
    return this;
  };

  /**
   * equals 
   * 
   * @param {simplexn.SimplicialComplex} simpcomp
   * @return {Boolean} true if this is equals to the given mplicial complex, false otherwise.
   * @api public
   */

  simplexn.SimplicialComplex.prototype.equals = function (simpcomp) {
    var pointset1 = this.pointset;
    var pointset2 = simpcomp.pointset;
    var topology1 = this.topology;
    var topology2 = simpcomp.topology;

    if (! pointset1.equals(pointset2)) return false;
    if (! topology1.equals(topology2)) return false;

    return true;
  };

  /**
   * extrude
   * 
   * @param {Array|Float32Array} hlist which must be made by positive numbers 
   *   or by an alternation of positive and negative numbers
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api private
   */

  simplexn.SimplicialComplex.prototype.extrude = function (hlist) {
    var hlist = hlist || [1];
    var hlistLength = hlist.length;
    var hlist0 = hlist[0];
    var hlist0isNegative = hlist0 < 0;
    var positiveQuotes = hlist.filter(function (h) {return h >= 0}).length;
    
    var oldRn = this.rn;
    var newRn = oldRn + 1;
    var pointset = this.pointset;
    var pointsetSize = pointset.size;
    var oldEmbeddedPoints;
    var oldPointsLength = pointsetSize * oldRn;
    var newPointsLength = pointsetSize * newRn;
    var newPoints = new Float32Array(newPointsLength);

    var newPointsetSize = (hlistLength + 1 - (hlist0isNegative)) * pointsetSize;
    var newPointset = new PointSet(newPointsetSize, newRn);

    var oldDim = this.dim;
    var newDim = oldDim + 1; 
    var topology = this.topology;
    var cellLength = oldDim + 1;
    var complex = topology.complexes[oldDim];
    var complexLength = complex.length;
    var complexSize = complexLength / cellLength;

    var newCellLength = newDim + 1;
    var newComplexLength = positiveQuotes * complexSize * newDim * newCellLength;
    var newComplex = new Uint32Array(newComplexLength);

    var tempLength = 2 * cellLength;
    var temp = new Uint32Array(tempLength);
    var tempIndx;
    var cIndx = 0 ;
    var exchange1, exchange2;
    var end;
    var quote = 0;
    var h, v, c, i, j;

    for (i = 0; i < complexLength; i += cellLength) {
      _quickSort(complex.subarray(i, i + cellLength));
    }

    this.embed();
    oldEmbeddedPoints = this.pointset.points;
    newPoints.set(oldEmbeddedPoints);
    if (!hlist0isNegative) newPointset.set(oldEmbeddedPoints);

    for (h = 0; h < hlistLength; h += 1) {
      quote += abs(hlist[h]);

      // add new points
      for (v = newRn - 1; v < newPointsLength; v += newRn) {
        newPoints[v] = quote;
      }
      newPointset.set(newPoints, (h + 1 - (hlist0isNegative)) * pointsetSize);

      // create new cells
      if (hlist[h] >= 0) {
        for (c = 0; c < complexSize; c += 1) {
          // fill temp with selected indexes
          for (i = 0; i < cellLength; i++) {
            tempIndx = complex[c*cellLength+i] + (h - (hlist0isNegative)) * pointsetSize;
            temp[i] = tempIndx;
            temp[i+cellLength] = tempIndx + pointsetSize;
          }
          
          // pick cells from temp, cellLength by cellLength
          for (i = 0; i < cellLength; i += 1) {
            end = i + cellLength + 1;
            for (j = i; j < end; j++) {
              newComplex[cIndx++] = temp[j];
            }
            // take care of orientation
            if (((newDim & 1) * (c) + (oldDim & 1) * i) & 1) {
              exchange1 = cIndx - 1;
              exchange2 = exchange1 - 1;
              _swap(newComplex, exchange1, exchange2);
            }
          }
        }
      }
    }

    this.pointset = newPointset;
    this.topology = new Topology(newComplex, newDim);
    return this;
  };

  /**
   * getFacet
   * Return the `dim` and `index` facet
   *
   * @param {Number} dim
   * @param {Number} index
   * @return {Array} the facet represented by a `dim+1` length array
   * @api private
   */

  simplexn.SimplicialComplex.prototype.getFacet = function (dim, index) {
    var topology = this.topology;
    var pointset = this.pointset;
    var cells = this.topology.complexes[dim];
    var size = dim + 1;
    var start = size*index;
    var cell = cells.subarray(start, start + size);
    var facet = [];
    var i;

    for (i = 0; i < size; i += 1) {
      facet.push(pointset.get(cell[i]));
    }

    return facet;
  };

  /**
   * centroids
   * Return a PointsSet of centroids of dim-cells.
   * 
   * @param {Number} [dim=this.dim]
   * @return {simplexn.PointSet}
   * @api public
   */

  simplexn.SimplicialComplex.prototype.centroids = function (dim) {
    var dim = dim || this.dim;
    var cellSize = dim + 1;
    var cellsSize = this.topology.complexes[dim].length / cellSize;
    var centroids = new PointSet(cellsSize, this.rn);
    var i;
    var centroid;

    for (i = 0; i < cellsSize; i += 1) {
      centroid = vectorAvg(this.getFacet(dim, i));
      centroids.set(centroid, i);
    }

    return centroids;

  };

  /**
   * explode
   *
   * @param {Array|Float32Array} values
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.explode = function (values) {
    var dim = this.dim;
    var values = values || [];
    var cell;
    var cellSize = dim + 1;
    var cells = this.topology.complexes[dim];
    var cellsSize = cells.length / cellSize;
    var newCell, newCells = [];
    var pointset = this.pointset;
    var newPointset = new PointSet(cellsSize*cellSize, this.rn);
    var centroids = this.centroids();
    var centroid;
    var translatedCentroid;
    var translationVect;
    var c, i;
    var indx = 0;

    for (c = 0; c < cellsSize; c += 1) {
      cell = cells.subarray(c*cellSize, c*cellSize + cellSize);
      newCell = [];
      centroid = centroids.get(c);
      translatedCentroid = vectorMul(centroid, values);
      translationVect = vectorSub(translatedCentroid, centroid);
      for (i = 0; i < cellSize; i += 1) {
        newCell.push(indx);
        newPointset.set(vectorAdd(pointset.get(cell[i]), translationVect), indx);
        indx++;
      }
      newCells.push(newCell);
    }

    this.pointset = newPointset;
    this.topology = new Topology(newCells);
    return this.merge();
  };

  /**
   * skeleton
   * 
   * @param {Number} dim
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.skeleton = function (dim) {
    this.topology.skeleton(dim);
    return this;
  };

  /**
   * boundary
   * 
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.boundary = function () {
    this.topology.boundary();
    return this;
  };

  /**
   * prod
   * Execute product of this simplicial complex 
   * by a given simplicial complex.
   * At the moment it's customed and tested only for following cases:
   * - 1-rn x 1-rn
   * - 1-rn x 2-rn
   * - 2-rn x 1-rn
   * 
   * @param {simplexn.SimplicialComplex} simpcomp
   * @return {simplexn.SimplicialComplex} this for chaining
   * @api public
   */

  simplexn.SimplicialComplex.prototype.prod = function(simpcomp) {
    if (this.rn > 1 && simpcomp.rn > 2) return;
    if (this.rn > 2 && simpcomp.rn > 1) return;

    var n = simpcomp.size - 1;
    var pointset = this.pointset.clone().prod(simpcomp.pointset);
    var quotes = [];

    while (n--) quotes.push(1);

    this.extrude(quotes);

    this.pointset = pointset;

    return this;
  };

  /**
   * Struct
   * A simplexn.Struct is a collection of simplexn.SimplicialComplex
   *
   * @constructor
   * @param {Array} items simplexn.SimplicialComplex or simplexn.Struct instaces;
   * @api public
   */

  var Struct = 
  simplexn.Struct = function (items) {
    var items = items || [];
    var complexes = [];
    var structs = [];

    items.forEach(function (item) {
      if (item instanceof SimplicialComplex) {
        complexes.push(item);
      } else if (item instanceof Struct) {
        structs.push(item);
      }
    });

    this.complexes = complexes;
    this.structs = structs;
  };

  /**
   * rotate
   * 
   * @param {Array|Uint32Array} dims
   * @param {Number} angle
   * @return {simplexn.Struct} this for chaining
   * @api public
   */

  simplexn.Struct.prototype.rotate = function (dims, angle) {
    this.complexes.forEach(function (complex) {
      complex.rotate(dims, angle);
    });

    this.structs.forEach(function (struct) {
      struct.rotate(dims, angle);
    });

    return this;
  };

  /**
   * scale
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {simplexn.Struct} this for chaining
   * @api public
   */

  simplexn.Struct.prototype.scale = function (dims, values) {
    this.complexes.forEach(function (complex) {
      complex.scale(dims, values);
    });

    this.structs.forEach(function (struct) {
      struct.scale(dims, values);
    });

    return this;
  };

  /**
   * translate
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {simplexn.Struct} this for chaining
   * @api public
   */

  simplexn.Struct.prototype.translate = function(dims, values) {
    this.complexes.forEach(function (complex) {
      complex.translate(dims, values);
    });

    this.structs.forEach(function (struct) {
      struct.translate(dims, values);
    });

    return this;
  };

  /**
   * clone
   * 
   * @return {simplexn.Struct} cloned Struct
   * @api public
   */

  simplexn.SimplicialComplex.prototype.clone = function () {
    var clone = new SimplicialComplex();
    clone.pointset = this.pointset.clone();
    clone.topology = this.topology.clone();
    return clone;
  };

  /**
   * geometries
   */

  simplexn.geometries = {};

  /**
   * simplex
   * 
   * @param {number} d
   * @return {simplexn.SimplicialComplex} a simplex
   * @api public
   */

  var simplex =
  simplexn.geometries.simplex = function (d) {
    var d = d !== undefined ? d + 1 : 1;
    var dim = d;
    var points0 = [];
    var points = [];
    var cells = [];

    while (d--) {
      points0.push(0);
      cells.unshift(d);
    }

    points =  matrixIdentity(dim);
    points.unshift(points0);

    return new SimplicialComplex(points, [cells]);
  }

  /**
   * polyline
   * 
   * @param {Array} points
   * @return {simplexn.SimplicialComplex} a simplex
   * @api public
   */

  var polyline =
  simplexn.geometries.polyline = function (points) {
    var points = points || [[]];
    var n = points.length - 1;
    var cells = [];
    var i;

    for (i = 0; i < n; i += 1) {
      cells.push([i, i+1]);
    }

    return (new SimplicialComplex(points, cells)).merge();
  };

  /**
   * polypoint
   * 
   * @param {Array} points
   * @return {simplexn.SimplicialComplex} a simplex
   * @api public
   */

  var polypoint =
  simplexn.geometries.polypoint = function (points) {
    var points = points || [[]];
    var cells = points.map(function (p,i) { 
      return [i];
    });

    return (new SimplicialComplex(points, cells));
  };

  /**
   * simplexGrid
   * 
   * @param {Array} quotesList is a list of hlist which must be made by positive numbers 
  *                 or by an alternation of positive and negative numbers
   * @return {simplexn.SimplicialComplex} a grid of simplexes
   * @api public
   */

  var simplexGrid = 
  simplexn.geometries.simplexGrid = function (quotesList) {
    var quotesList = quotesList ? quotesList.slice(0) : [[]];
    var quotesListHead = quotesList.shift();
    var quotesListHead0 = quotesListHead[0];
    var quotesListHead0isNeg = quotesListHead0 <= 0;
    var points = quotesListHead0isNeg ? [] : [[0]];
    var length = quotesList.length;
    var complex = [];
    var simpcomp;
    var quote = 0;
    var indx;

    quotesListHead.forEach(function (height, i) {
      quote += abs(height);
      points.push([quote]);
      if (height > 0) {
        indx = i - (quotesListHead0isNeg);
        complex.push([indx,indx+1]);
      }
    });

    simpcomp = new SimplicialComplex(points, complex);

    quotesList.forEach(function (quotes) {
      simpcomp.extrude(quotes);
    });
  
    return simpcomp;
  };

  /**
   * cuboid
   * 
   * @param {Array} sideds
   * @return {simplexn.SimplicialComplex} a cuboidal simplicial complex
   * @api public
   */

  simplexn.geometries.cuboid = function (sides) {
    sides = sides.map(function (s) { return [s]; });
    return simplexGrid(sides);
  };

  /**
   * intervals
   *
   * @param {Array} values
   * @return {simplexn.SimplicialComplex} intervals
   * @api public
   */

  var intervals = 
  simplexn.geometries.intervals = function (tip, n) {
    var values = [];
    var value = tip/n;
    
    while (n--) values.push(value);

    return simplexGrid([values]);
  };

  /**
   * domain
   *
   * @param {Array} ends
   * @param {Array} ns
   * @return {simplexn.SimplicialComplex} domain
   * @api public
   */

  var domain = 
  simplexn.geometries.domain = function (ends, ns) {
    var ends = ends || [0, 2*Math.PI];
    var ns = ns || [36];
    var length = ends.length;
    var endsn = ends[0];
    var begin = endsn[0];
    var end = endsn[1];
    var domain = intervals(end - begin, ns[0]).translate([0], [begin]);
    var values;
    var value;
    var i;
    var n;

    for (i = 1; i < length; i += 1) {
      endsn = ends[i];
      begin = endsn[0];
      end = endsn[1];
      n = ns[i];
      values = [];
      value = (end - begin)/n;

      while (n--) values.push(value);

      domain = domain.extrude(values).translate([i], [begin]);
    }

    return domain;
  };

  /**
   * cube
   * 
   * @param {Number} d
   * @return {simplexn.SimplicialComplex} a dim-dimendional cube
   * @api public
   */

  simplexn.geometries.cube = function (d) {
    var d = d || 1;
    var quotes = [];
    while (d--) quotes.push([1]);
    return simplexGrid(quotes);
  };

  /**
   * circle
   * 
   * @param {Number} [radius=1]
   * @param {Number} [n=32] 
   * @return {simplexn.SimplicialComplex} a circle
   * @api public
   */

  simplexn.geometries.circle = function (radius, n) {
    var r = radius || 1;
    var n = n || 32;
    var domain = intervals(2 * pi, n);
    
    domain.map(function (v) { 
      return [r * sin(v[0]), r * cos(v[0])]; 
    }, true);

    return domain;
  };

  /**
   * disk
   * 
   * @param {Number} [radius=1]
   * @param {Number} [n=32]
   * @param {Number} [m=1] 
   * @return {simplexn.SimplicialComplex} a disk
   * @api public
   */

  simplexn.geometries.disk = function (radius, n, m) {
    var radius = radius || 1;
    var n = n || 32;
    var m = m || 1;
    var nQuote = 2 * pi / n;
    var mQuote = radius / m;
    var nQuotes = [];
    var mQuotes = [];
    var domain;

    while (n--) nQuotes.push(nQuote);
    while (m--) mQuotes.push(mQuote);
    domain = simplexGrid([nQuotes,mQuotes]);
    domain.map(function (coords) {
      var u = coords[0];
      var v = coords[1];
      return [v*sin(u), v*cos(u)];
    }, true);

    return domain;
  };

  /**
   * cylinderSurface
   * Produces a cylindrical surface of radius r and heigth h.
   * 
   * @param {Number} [r=1]
   * @param {Number} [h=1]
   * @param {Number} [n=16]
   * @param {Number} [m=2] 
   * @return {simplexn.SimplicialComplex} a cylindrical surface
   * @api public
   */
  simplexn.geometries.cylinderSurface = function (r, h, n, m) {
    var r = r || 1;
    var h = h || 1;
    var n = n || 16;
    var m = m || 2;
    var domain = simplexGrid([_repeat(2*pi/n, n), _repeat(1./m, m)]);

    domain.map(function(v) {
      return [
        r * cos(v[0])
      , r * sin(v[0])
      , h * v[1]
      ];
    }, true);

    return domain;
  };

  /**
   * cylinderSolid
   * Produces a solid cylindrer with radius r and heigth h.
   *
   * @param {Number} [R=1] 
   * @param {Number} [r=0]
   * @param {Number} [h=1]
   * @param {Number} [n=16]
   * @param {Number} [m=1]
   * @param {Number} [p=1] 
   * @return {simplexn.SimplicialComplex} a cylinder
   * @api public
   */

  var cilinderSolid =
  simplexn.geometries.cylinderSolid = function (R, r, h, n, m, p) {
    var R = R || 1.; 
    var r = r || 0.; 
    var h = h || 1.;
    var n = n || 16;
    var m = m || 1;
    var p = p || 1; 
    var domain = simplexGrid([_repeat(2*pi/n, n), _repeat((R-r)/m, m), _repeat(h/p, p)]);
    
    domain.translate([1],[r]).map(function(v) {
      return [
        v[1] * sin(v[0])
      , v[1] * cos(v[0])
      , v[2]
      ];
    }, true);

    return domain;
  };

  /**
   * torusSurface
   *
   * produces a toroidal surface of radiuses r,R 
   * approximated with n x m x 2 triangles
   *
   * @param {Number} [r=1] r
   * @param {Number} [R=3] R
   * @param {Number} [n=12] n
   * @param {Number} [m=8] m
   * @return {simplexn.SimplicialComplex} torus surface
   * @api public
   */

  var torusSurface =
  simplexn.geometries.torusSurface = function (r, R, n, m) {
    var r = r || 0.5;
    var R = R || 1.5;
    var n = n || 12;
    var m = m || 8;
    var domain = simplexGrid([ _repeat(2*pi/n, n), _repeat(2*pi/m, m)]);

    domain.map(function (v) {
      return [
          (R + r * cos(v[1])) * cos(v[0])
        , (R + r * cos(v[1])) * sin(v[0])
        , r * sin(v[1])
      ];
    }, true);

    return domain;
  };

 /**
  * torusSolid
  *
  * produces a toroidal surface of radiuses r,R 
  * approximated with n x m x 2 triangles
  *
  * @param {Number} [r=1] r
  * @param {Number} [R=3] r
  * @param {Number} [n=12] n
  * @param {Number} [m=8] m
  * @param {Number} [p=8] p
  * @return {simplexn.SimplicialComplex} torus solid
  * @api public
  */

  var torusSolid = 
  simplexn.geometries.torusSolid = function (r, R, n, m, p) {
    var r = r || 1;
    var R = R || 3;
    var n = n || 12;
    var m = m || 8;
    var p = p || 2;
    var domain = simplexGrid([ _repeat(2*pi/n, n), _repeat(2*pi/m, m), _repeat(1/pi, p)]);

    domain.map(function (v) {
      return [
          (R + r * v[2] * cos(v[0])) * cos(v[1])
        , (R + r * v[2] * cos(v[0])) * -sin(v[1])
        , r * v[2] * sin(v[0])
      ];
    }, true);

    return domain;
  };

  /**
   * triangleStrip
   * 
   * @param {Array} points
   * @return {simplexn.SimplicialComplex} triangle strip
   * @api public
   */

  var triangleStrip = 
  simplexn.geometries.triangleStrip = function (points) {
    var n = points.length;
    var cells = [];
    var i;
    
    for (i = 2; i < n; i += 1) {
      if (cells.length & 1) {
        cells.push([i-1, i-2, i-0]);
      }
      else {
        cells.push([i-2, i-1, i-0]);
      }
    }

    return new SimplicialComplex(points, cells);
  };


  /**
   * triangleFan
   * 
   * @param {Array} points
   * @return {simplexn.SimplicialComplex} triangle strip
   * @api public
   */

  var triangleFan = 
  simplexn.geometries.triangleFan = function (points) {
    var n = points.length;
    var cells = [];
    var i;
    
    for (i = 2; i < n; i += 1) {
      cells.push([0, i-1, i]);
    }

    return new SimplicialComplex(points, cells);
  };

  /**
   * helix
   *
   * @param {Number} [r=1] r
   * @param {Number} [pitch=1] pitch
   * @param {Number} [n=24] n
   * @param {Number} [turns=1] turns
   * @return {simplexn.SimplicialComplex} helix
   * @api public
   */

  var helix = 
  simplexn.geometries.helix = function (r, pitch, n, turns) {
    var r = r || 1;
    var pitch = pitch || 1;
    var n = n || 24;
    var turns = turns || 8;
    var domain = intervals(2*pi*turns, n*turns);

    domain.map(function (v) {
      return [
          r * sin(v[0])
        , r * cos(v[0])
        , pitch / (2*pi) * v[0]
      ];
    }, true);

    return domain;
  };

}(this));/*!
 * ƒ
 * JavaScript functional library
 * Copyright (c) 2012 Enrico Marino <enrico.marino@email.com> (onirame.no.de)
 * MIT License
 */

 !(function (exports) {

  /**
   * Library namespace.
   */

   var ƒ = exports.ƒ = exports.f = {};

  /**
   * Library version.
   */

  ƒ.version = '0.6.5';

  /**
   * choose 
   * binomial coefficients
   * 
   * @example 
   *   choose([7,5]); //21
   *
   * @param {Array} pair
   * @return {Function} binomial coefficient
   * @api public
   */

  ƒ.choose = function (pair) {
    var n = pair[0];
    var k = pair[1];
    var coeff = 1;
    var i;

    for (i = n-k+1; i <= n; i++) {
      coeff *= i;
    }

    for (i = 1; i <= k; i++) {
      coeff /= i;
    }

    return coeff;
  };

  /**
   * sel 
   * select n-th element from array
   * 
   * @example 
   *   sel(2)([0,1,2,3]); //2
   *
   * @param {Number} index of selection
   * @return {Function}
   *   @param {Array} array
   *   @return {*} array[index]
   * @api public
   */
  
  ƒ.sel = function (index) {
    return function (array) {
      return array[index];
    };
  };

  ƒ.s0 = ƒ.sel(0);
  ƒ.s1 = ƒ.sel(1);
  ƒ.s2 = ƒ.sel(2);
  ƒ.s3 = ƒ.sel(3);
  ƒ.s4 = ƒ.sel(4);
  ƒ.s5 = ƒ.sel(5);
  ƒ.s6 = ƒ.sel(6);
  ƒ.s7 = ƒ.sel(7);
  ƒ.s8 = ƒ.sel(8);
  ƒ.s9 = ƒ.sel(9);

  /**
   * apply 
   * apply([f,x]) 
   * apply `f` to `x`.
   * 
   * @example 
   *   apply([Math.cos, Math.PI/3]); //0.5
   *
   * @param {Function} [f = pair[0]] the function to apply
   * @param {*} [x = pair[1]] the value to apply `f` to.
   * @return {*} the result of `f(x)`.
   * @api public
   */
  
  ƒ.apply = function (pair) {
    var f = pair[0];
    var x = pair[1];
    return f(x);
  };

  /**
   * aa
   * aa(f)(array)
   * apply `f` to each element of `array`.
   * 
   * @example aa(function (x) {return x * 2;})([1,3,5,7,9]); //[2,6,10,14,18]
   *
   * @param {Function} function f
   * @return {Function}
   *    @param {Array} array [a0,a1,...,an]
   *    @return {Array} [f(a0),f(a1),...,f(an)]
   * @api public
   */

  ƒ.aa = function (f) {
    return function (array) {
      return array.map(function (element) {
        return f(element);
      });
    };
  };

  /**
   * comp2
   * returns the composition of the given functions
   * 
   * @example
   *   comp2([
   *     function (x) {return x * 2;}, 
   *     function (y) {return y - 1;}
   *   ])(5); //8
   * 
   * @param {Array} functions array of functions to compose
   * @param {Function} [functions[0]] f
   * @param {Function} [functions[1]] g
   * @return {Function} the composition of the given functions
   * @api public
   */

  ƒ.comp2 = function (functions) {
    var f = functions[0];
    var g = functions[1];
    return function (x) {
      return f(g(x));
    };
  };

  /**
   * comp
   * returns the composition of the given functions
   * 
   * @example
   *   comp([
   *     function (x) {return x + 1;}, 
   *     function (y) {return y * 2;},
   *     function (z) {return z - 1;}
   *   ])(5); //3
   *
   * @param {Array} functions array of functions to compose
   * @return {Function} the composition of the given functions
   * @api public
   */

  ƒ.comp = function (functions) {
    return functions.reduce(function (f, g) {
      return function (x) {
        return f(g(x));
      };
    });
  };

  /**
   * cons
   * Apply each function of the given array `functions` to the given value `x`,
   * and return the array of application values
   *
   * @example
   *   cons([
   *     function (x) {return x - 1;}, 
   *     function (y) {return y * 2;},
   *     function (z) {return z % 3}
   *   ])(5); //[4,10,2]
   *
   * @param {Array} functions
   * @return {Array} the array of application values
   * @api public
   */

  ƒ.cons = function (functions) {
    return function (x) {
      return functions.map(function (f) {
        return f(x);
      });
    };
  };

  /**
   * id
   * return the given `value`
   * 
   * @param value
   * @return the given `value`
   * @api public
   */

  ƒ.id = function (value) {
    return value;
  };
  
  /**
   * k
   * return a function that return the given `value`
   * 
   * @param value
   * @return {Function} a function that return `value`
   * @api public
   */

  ƒ.k = function (value) {
    return function () {
      return value;
    };
  };

  /**
   * cat
   * catenates `args`, an array of arrays, by eliminating a level of nesting
   *
   * @example
   *   cat([
   *     [0,1,2],
   *     [3,4,5,6],
   *     [7,8,9,10,11]
   *   ]); //[0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11]
   *
   * @param {Array} arrays array of arrays
   * @return {Array} array created eliminating a level of nesting
   * @api public 
   */
  
  ƒ.cat = function (arrays) {
    var result = [];
    arrays.forEach(function (array) {
      result = result.concat(array);
    });
    return result;
  };

  /**
   * distl
   * distribute left: 
   * returns the `pair` sequence with `value` and the elements of `array` 
   * 
   * @example
   *   distl(['a',[0,1,2,3,4]]); //[['a',0],['a',1],['a',2],['a',3],['a',4]]
   *
   * @param {Array} pair
   * @param {Array} [pair[0]] array
   * @param {Any}   [pair[1]] value
   * @return the `pair` sequence with `value` and the elements of `array`
   * @api public
   */

  ƒ.distl = function (pair) {
    var value = pair[0];
    var array = pair[1];
    return array.map(function (item) {
      return [value, item];
    });
  };

  /**
   * distr
   * distribute right: 
   * returns the `pair` sequence with the elements of `array` and `value`
   * 
   * @example
   *   distr([[0,1,2,3,4],'a']); //[[0,'a'],[1,'a'],[2,'a'],[3,'a'],[4,'a']]
   *
   * @param {Array} pair
   * @param {Array} [pair[0]] array
   * @param {Any}   [pair[1]] value
   * @return the `pair` sequence with the elements of `array` and `value`
   * @api public
   */

  ƒ.distr = function (pair) {
    var array = pair[0];
    var value = pair[1];
    return array.map(function (item) {
      return [item, value];
    });
  };

  /**
   * insl
   * insert left operator  
   * given a binary associative `operator` 
   * returns a function that given an array 
   * returns the riduction of the array by the operator.
   * 
   * @param {Function} operator binary operator
   * @return {Function} function that apply `operator` to the given `array`
   * @api public
   */

  ƒ.insl = function (operator) {
    return function (array) {
      return array.reduce(operator);
    };
  };

  /**
   * insr
   * insert right operator  
   * given a binary associative `operator` 
   * returns a function that given an array 
   * returns the right riduction of the array by the operator.
   * 
   * @param {Function} operator binary operator
   * @return {Function} function that apply `operator` to the given `array`
   * @api public
   */

  ƒ.insr = function (operator) {
    return function (array) {
      return array.reduceRight(operator);
    };
  };
  
  /**
   * al
   * append left
   * append `item` on the left of `array`
   * 
   * @param {Array} pair
   * @param {Array} [pair[0]] item
   * @param {Any}   [pair[1]] array
   * @return {Array} `array` concatenated with `item`
   * @api public
   */

  ƒ.al = function (pair) {
    var item = pair[0];
    var array = pair[1];
    return [item].concat(array);
  };

  /**
   * ar
   * append right
   * append `item` on the right of `array`
   * 
   * @param {Array} pair
   * @param {Any}   [pair[0]] array
   * @param {Array} [pair[1]] item
   * @return {Array} `array` concatenated with `item`
   * @api public
   */

  ƒ.ar = function (pair) {
    var array = pair[0];
    var item = pair[1];
    return array.concat([item]);
  };

  /**
   * last
   * returns the last element of the given `array`
   *
   * @example
   *   last([0,1,2,3,4,5]); //5
   *
   * @param {Array} array array
   * @return {Any} the last element of `array`
   * @api public
   */

  ƒ.last = function (array) {
    return array[array.length - 1];
  };

  /** 
   * list
   * returns an array containing `arg`. 
   *
   * @param {Array} arg
   * @return {Array} array containing `arg`
   * @api public
   */

  ƒ.list = function (arg) {
    return [arg];
  };

  /**
   * len
   * returns the length of the given `array`
   * 
   * @param {Array} array array
   * @return {Number} the length of the given `array`
   * @api public
   */

  ƒ.len = function (array) {
    return array.length;
  };

  /**
   * first
   * returns the first element of the given `array`.
   * 
   * @param {Array} array array
   * @return the first element of the given `array`.
   * @api public
   */ 

  ƒ.first = function (array) {
    return array[0];
  };

  /** 
   * reverse
   * returns the given `array` in reverse order
   *
   * @example
   *   reverse([0,1,2,3,4,5]); //[5,4,3,2,1,0]
   *
   * @param {Array} array array
   * @return {Array} the given `array` in reverse order
   * @api public
   */

  ƒ.reverse = function (array) {
    var result = [];
    var i;
    for (i = array.length - 1; i >= 0; i--) {
      result.push(array[i]);
    }
    return result;
  };
  
  /**
   * tail
   * returns the non-empty `array` but its `first` element
   * 
   * @example
   *   tail([0,1,2,3,4,5]); //[1,2,3,4,5]
   *
   * @param {Array} array array
   * @return {Array} the tail of the given `array`
   * @api public
   */

  ƒ.tail = function (array) {
    return array.slice(1);
  };

  /**
   * butlast
   * returns the non-empty `array` but its `last` element
   * 
   * @example
   *   butlast([0,1,2,3,4,5]); //[0,1,2,3,4]
   * @example
  *    butlast([]); //[]
   *
   * @return {Array} the non-empty `array` but its `last` element
   * @api public
   */

  ƒ.butlast = function (array) {
    return array.slice(0,-1);
  };

  /**
   * repeat
   * returns an array with `n` repetitions of `value`
   * 
   * @example
   *   repeat(3)(12); //[12,12,12]
   *
   * @param {Number} n number of repetitions
   * @return {Function} a function that given `value` 
   *   returns an array with `n` repetitions of `value`
   * @api public
   */

  ƒ.repeat = function (n) {
    return function (value) {
      var result = [];
      var i;
      for (i = 0; i < n; i++) {
        result.push(value);
      }
      return result;
    };
  };

  /**
   * replica
   * repeat list and catenate.
   * 
   * @example
   *   replica(3)(['A',1]); //["A", 1, "A", 1, "A", 1]
   *
   * @param {Number} n number of repetitions
   * @return {Function} a function that given `value` 
   *   returns an array with `n` repetitions of `value` concatenated
   * @api public
   */
  
  ƒ.replica = function (n) {
    return function (value) { 
      var result = [];
      var i;
      for (i = 0; i < n; i++) {
        result = result.concat(value);
      }
      return result;
    };
  };

  /**
   * bigger
   * binary operator that returns the greater of the given pair
   * 
   * @example
   *   bigger([4,9]); //9
   *
   * @param {Array} pair
   * @return {Number} return the greater of the `pair`
   * @api public
   */

  ƒ.bigger = function (pair) {
    var a = pair[0];
    var b = pair[1];
    return a > b ? a : b;
  };

  /**
   * smaller
   * binary operator that returns the smaller of the given pair
   * 
   * @example
   *   bigger([4,9]); //9
   *
   * @param {Array} pair
   * @return {Number} return the smaller of the `pair`
   * @api public
   */

  ƒ.smaller = function (pair) {
    var a = pair[0];
    var b = pair[1];
    return a < b ? a : b;
  };

  /**
   * biggest
   * returns the greatest of the given `values`
   * 
   * @example
   *   biggest([4,9,2,8,1,7]); //9
   *
   * @param {Array} values values
   * @return {Number} return the greatest of the given `values`
   * @api public
   */

  ƒ.biggest = function (values) {
    return Math.max.apply(null, values);
  };

  /**
   * smallest
   * returns the smallest of the given `values`
   *
   * @example
   *   smallest([4,9,2,8,1,7]); //1
   *
   * @param {Array} values values
   * @return {Number} return the smallest of the given `values`
   * @api public
   */

  ƒ.smallest = function (values) {
    return Math.min.apply(null, values);
  };

  /**
   * sum
   * returns the sum of the given `values` (numbers or arrays)
   * 
   * @example
   *   sum([[1,2,3],[2,3,4],[3,4,5]]); //[6,9,12]
   * @example
   *   sum([1,2,3,4]); //10
   *
   * @param {Array} values values
   * @return {Number} return the sum of the given `values`
   * @api public
   */

  ƒ.sum = function (values) {
    if (values[0] instanceof Array) {
      return values.reduce(function (prev, curr) {
        return prev.map(function (value, i) {
          return value + curr[i];
        });
      });
    }
      
    return values.reduce(function (prev, curr) {
      return prev + curr;
    });
  };

  /**
   * sub
   * returns the difference of the given `values` (numbers or arrays)
   * 
   * @param {Array} values values
   * @return {Number} return the difference of the given `values`
   * @api public
   */

  ƒ.sub = function (values) {
    if (values[0] instanceof Array) {
      return values.reduce(function (prev, curr) {
        return prev.map(function (value, i) {
          return value - curr[i];
        });
      });
    }
      
    return values.reduce(function (prev, curr) {
      return prev - curr;
    });
  };
  
  /**
   * mul
   * multiplicate the given `values`
   *
   * @example
   *   mul([[1,2,3],[2,3,4],[3,4,5]]); //[6,24,60]
   * @example
   *   mul([1,2,3,4]); //24
   *
   * @param {Array} values values
   * @return {Number} return the product of the given `values`
   * @api public
   */

  ƒ.mul = function (values) {
    if (values[0] instanceof Array) {
      return values.reduce(function (prev, curr) {
        return prev.map(function (value, i) {
          return value * curr[i];
        });
      });
    }
      
    return values.reduce(function (prev, curr) {
      return prev * curr;
    });
  };

  /**
   * div
   * divide the given `values`
   *    
   * @param {Array} values values
   * @return {Number} return the division of the given `values`
   * @api public
   */

  ƒ.div = function (values) {
    if (values[0] instanceof Array) {
      return values.reduce(function (prev, curr) {
        return prev.map(function (value, i) {
          return value / curr[i];
        });
      });
    }
      
    return values.reduce(function (prev, curr) {
      return prev / curr;
    });
  };
  
  /** 
   * cart
   * cartesian product
   * 
   * @example
   *  cart([[1,2],['a','b']]); //[[1,"a"],[1,"b"],[2,"a"],[2,"b"]]
   *
   * @param {Array} args
   * @return {Array} the cartesian product of `args`
   * @api public
   */

  ƒ.cart = function (args) {
    return args.reduce(function (a, b) {
      var ret = [];
      a.forEach(function (a) {
        b.forEach(function (b) {
          ret.push(a.concat([b]));
        });
      });
      return ret;
    }, [[]]);
  };

  /**
   * Math constant
   */

  /**
   * E
   *
   * @constant
   */

  ƒ.E = Math.E;

  /**
   * LN2
   *
   * @constant
   */

  ƒ.LN2 = Math.LN2;

  /**
   * LN10
   *
   * @constant
   */

  ƒ.LN10 = Math.LN10;

  /**
   * LOG2E
   *
   * @constant
   */

  ƒ.LOG2E = Math.LOG2E;

  /**
   * LOG10E
   *
   * @constant
   */

  ƒ.LOG10E = Math.LOG10E;

  /**
   * PI
   *
   * @constant
   */

  ƒ.PI = Math.PI;

  /**
   * SQRT1_2
   *
   * @constant
   */

  ƒ.SQRT1_2 = Math.SQRT1_2;

  /**
   * SQRT2
   *
   * @constant
   */

  ƒ.SQRT2 = Math.SQRT2;

  /**
   * abs
   * returns the absolute value of the given number
   *
   * @param {Number} number
   * @return {Number} the absolute value of the given number
   * @api public
   */

  ƒ.abs = Math.abs;

  /**
   * acos
   * returns the arc-cosine of the given number
   *
   * @param {Number} number
   * @return {Number} the arc-cosine of the given number
   * @api public
   */

  ƒ.acos = Math.acos;

  /**
   * asin
   * returns the arc-sine of the given number
   *
   * @param {Number} number
   * @return {Number} the arc-sine of the given number
   * @api public
   */

  ƒ.asin = Math.asin;

  /**
   * atan
   * returns the arc-tangent of the given number
   *
   * @param {Number} number
   * @return {Number} the arc-tangent of the given number
   * @api public
   */

  ƒ.atan = Math.atan;

  /**
   * atan2
   * returns the squared arc-tangent of the given number
   *
   * @param {Number} number
   * @return {Number} the squred arc-tangent of the given number
   * @api public
   */

  ƒ.atan2 = Math.atan2; 

  /**
   * ceil
   * returns the ceil of the given number
   *
   * @param {Number} number
   * @return {Number} the ceil of the given number
   * @api public
   */

  ƒ.ceil = Math.ceil; 

  /**
   * cos
   * returns the cosine of the given number
   *
   * @param {Number} number
   * @return {Number} the cosine of the given number
   * @api public
   */

  ƒ.cos = Math.cos; 

  /**
   * exp
   * returns the exponential of the given value (e^value)
   *
   * @param {Array} pair
   * @param {Number} [pair[0]] value
   * @param {Number} [pari[1]] n
   * @return {Number} the `n`-th power of the given number
   * @api public
   */

  ƒ.exp = Math.exp;

  /**
   * floor
   * returns the floor of the given number
   *
   * @param {Number} value
   * @return {Number} the floor of the given number
   * @api public
   */

  ƒ.floor = Math.floor;

  /**
   * log
   * returns the log of the given number
   *
   * @param {Number} value
   * @return {Number} the log of the given number
   * @api public
   */

  ƒ.log = Math.log;

  /**
   * floor
   * returns the floor of the given number
   *
   * @param {Number} value
   * @return {Number} the floor of the given number
   * @api public
   */

  ƒ.floor = Math.floor;

  /**
   * power
   * returns the n-th power of the given value (value^n)
   *
   * @param {Array} pair
   * @param {Number} [pair[0]] value
   * @param {Number} [pari[1]] n
   * @return {Number} the `n`-th power of the given number
   * @api public
   */

  ƒ.pow = function (pair) {
    return Math.pow.apply(null, pair);
  }; 

  /**
   * random
   * returns a random number in [0, 1) interval
   *
   * @return {Number} a random number in [0, 1) interval
   * @api public
   */

  ƒ.random = Math.random;  

  /**
   * round
   * returns the given number rounded 
   *
   * @return {Number} the given number rounded
   * @api public
   */

  ƒ.round = Math.round;  

  /**
   * sin
   * returns the sine of the given number
   *
   * @return {Number} the sin of the given number
   * @api public
   */

  ƒ.sin = Math.sin;  

  /**
   * sqrt
   * returns the squared root of the given number
   *
   * @return {Number} the squared root of the given number
   * @api public
   */

  var sqrt = 
  ƒ.sqrt = Math.sqrt; 

  /**
   * tan
   * returns the tan of the given number
   *
   * @return {Number} the tan of the given number
   * @api public
   */

  ƒ.tan = Math.tan;

  /**
   * trans
   * transpose the given matrix
   *
   * @example
   *   trans([[0,1,2],[3,4,5],[6,7,8]]); //[[0,3,6],[1,4,7],[2,5,8]
   *
   * @param {Array} matrix matrix
   * @return {Array} the transpose of the given matrix
   * @api public
   */

  ƒ.trans = function (matrix) {
    var result = [];

    matrix.forEach(function (row, i) {
      row.forEach(function (value, j) {
        (result[j] = result[j] || [])[i] = value;
      });
    });

    return result;
  };

  /**
   * prod
   * product scalar by vector
   *
   * @example
   *   prod([2, [0,1,1]]); //[0,2,2]
   *
   * @param {Array} pair
   * @param {Number} [pair[0]] scalar
   * @param {Array} [pair[1]] vector
   * @return {Array} the given `vector` scaled by the given `scalar`
   * @api public
   */

  ƒ.prod = function (pair) {
    var scalar = pair[0];
    var vector = pair[1];
    var result = vector.map(function (value) {
      return scalar * value;
    });
    return result;
  };

  /**
   * innerprod
   * inner (or scalar) product
   *
   * @example
   *   innerprod([[0,1,1],[2,3,1]]); //4
   *
   * @param {Array} pair
   * @param {Array} [pair[0]] v1
   * @param {Array} [pair[1]] v2
   * @return {Number} the inner producto of the given vectors
   * @api public
   */

  ƒ.innerprod = function (pair) {
    var v1 = pair[0];
    var v2 = pair[1];
    var result = 0;
    v1.forEach(function (value, i) {
      result += value * v2[i];
    });
    return result;
  };

  /**
   * vectnorm
   * returns the norm of the given `vector`
   *
   * @example
   *   vectrnorm([0,3,4]); //5
   *
   * @param {Array} vector
   * @return {Array} the the norm of the given `vector`
   * @api public
   */

  var vectnorm = 
  ƒ.vectnorm = function (vector) {
    return sqrt(vector.reduce(function (prev, current) {
      return prev + current * current;
    }, 0));
  };

  /**
   * unitvect
   * returns the unit vector of the given `vector`
   *
   * @example
   *   univect([0,3,4]); //[0, 0.6, 0.8]
   *
   * @param {Array} vector
   * @return {Array} the the unit vector of the given `vector`
   * @api public
   */

  var unitvect = 
  ƒ.unitvect = function (vector) {
    var norm = vectnorm(vector);
    var result = vector.map(function (value) {
      return value / norm;
    });
    return result;
  };

  /**
   * matsum
   * matrix sum
   *
   * @example
   *   matsum([[[1,0],[0,1]],[[1,0],[[1,0]]]); //[[2,0],[1,1]]
   *
   * @param {Array} pair
   * @param {Array} [pair[0]] m1
   * @param {Array} [pair[1]] m2
   * @return {Array} the sum of the given matrices
   * @api public
   */

  var matsum = 
  ƒ.matsum = function (pair) {
    var m1 = pair[0];
    var m2 = pair[1];
    var result = [];
    m1.forEach(function (row, i) {
      result[i] = [];
      row.forEach(function (value, j) {
        result[i][j] = value + m2[i][j];
      });
    })
    return result;
  };

  /**
   * matprod
   * matrix product
   *
   * @param {Array} pair
   * @param {Array} [pair[0]] m1
   * @param {Array} [pair[1]] m2
   * @return {Array} the matrix product of of the given matrices
   * @api public
   */

  var matprod = 
  ƒ.matprod = function (pair) {
    var m1 = pair[0];
    var m2 = pair[1];
    var n = m1.length;
    var m = m1[0].length;
    var i;
    var j;
    var k;
    var result = [];

    for (i = 0; i < n; i += 1) {
      result[i] = [];
      for (j = 0; j < m; j += 1) {
        result[i][j] = 0;
      }
    }

    for (i = 0; i < n; i += 1) {
      for (j = 0; j < m; j += 1) {
        for (k = 0; k < m; k += 1) {
          result[i][j] += m1[i][k] * m2[k][j];
        }
      }
    }

    return result;
  };

  /**
   * identity
   * matrix identity
   *
   * @example
   *   identity(4); //[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
   *
   * @param {Number} n dim of the identity matrix
   * @return {Array} an identity n x n matrix
   * @api public
   */

  var identity = 
  ƒ.identity = function (n) {
    var result = [];
    var i;

    for (i = 0; i < n; i += 1) {
      result[i] = [];
      for (j = 0; j < n; j += 1) {
        result[i][j] = (i === j) ? 1 : 0;
      }
    }

    return result;
  };

  /**
   * vectprod
   * vector product
   *
   * @param {Array} pair
   * @param {Array} [pair[0]] u
   * @param {Array} [pair[1]] v
   * @return {Array} the vector product of of the given vectors
   * @api public
   */

  var vectprod = 
  ƒ.vectprod = function (pair) {
    var u = pair[0];
    var v = pair[1];
    var result = [];

    result[0] = u[1]*v[2] - u[2]*v[1];
    result[1] = u[2]*v[0] - u[0]*v[2];
    result[2] = u[0]*v[1] - u[1]*v[0];

    return result;
  };

  /**
   * mat
   * create a nxm matrix from the given vector
   *
   * @example
   *   mat([3,2])([0,1,2,3,4,5]); //[[0,1],[2,3],[4,5]]
   *
   * @param {Array} dims
   * @param {Array} [dims[0]] n number of rows
   * @param {Array} [pair[1]] m number of cols
   * @return {Function} 
   *   @param {Array} vector
   *   @return {Array} matrix nxm from `vector` values
   * @api public
   */

  var mat = 
  ƒ.mat = function (dims) {
    var n = dims[0];
    var m = dims[1];

    return function (vector) {
      var result = [];
      var i;

      for (i = 0; i < n; i += 1) {
        result[i] = vector.slice(i*m, i*m+m);
      }
      
      return result;
    };
  };

}(this));
/*!
 * plasm.js
 * JavaScript Programming Language for Solid Modeling
 * Copyright (c) 2012 cvd-lab <cvdlab@email.com> (https://github.com/cvd-lab/)
 * MIT License
 */

 !(function (exports) {

  var toString = {}.toString;
  var max = Math.max;
  var min = Math.min;
  var abs = Math.abs;
  var abs = Math.abs;

  /**
   * Library namespace.
   */

  var plasm = exports.plasm = {};

  /**
   * Library version.
   */

  plasm.version = '0.1.5';

  /*
   * @intro Create the plasm viewer.
   * @param {Element} container
   * @param {Element} inspector
   * @api public
   */

  var Plasm = exports.Plasm =
  plasm.Viewer = function (container, inspector) {
    if (!(this instanceof plasm.Viewer)) {
      return new plasm.Viewer(container);
    }
    // if (!Detector.webgl) {
    //   Detector.addGetWebGLMessage();
    // }
    if (typeof container === 'string') {
      container = document.getElementById(container);
    }
    if (typeof inspector === 'string') {
      inspector = document.getElementById(inspector);
    }

    var scene = this.scene = new plasm.Scene();

    var camera = this.camera = new plasm.Camera();
    scene.add(camera);

    var controls = this.controls = new plasm.Controls(camera, scene);

    var light = this.light = new THREE.AmbientLight(0xeeeeee);
    scene.root.add(light);

    var axes = this.axes();
    axes.draw();

    var engine = Detector.webgl ? THREE.WebGLRenderer : THREE.CanvasRenderer;
    var renderer = this.renderer = new engine({ antialias: true });
    renderer.setClearColorHex(0xefefef, 1);
    renderer.setSize(window.innerWidth, window.innerHeight);

    this.container = container;
    this.container.appendChild(this.renderer.domElement);

    if (inspector) {
      var stats = this.stats = new Stats();
      inspector.appendChild(stats.domElement);
    }

    function animate () {
      requestAnimationFrame(animate);
      controls.update();
      renderer.render(scene.root, camera.optics);
      if (inspector) stats.update();
    }

    function resize () {
      camera.optics.aspect = window.innerWidth / window.innerHeight;
      camera.optics.updateProjectionMatrix();
      renderer.setSize(window.innerWidth, window.innerHeight);
    }

    window.addEventListener('resize', resize, false);

    console.log('Plasm.js');

    animate();
  };

  /**
   * Create a scene.
   *
   * @return {plasm.Scene} scene
   * @api public
   */

  plasm.Scene = function () {
    this.root = new THREE.Scene();
  };

  /**
   * Add object to the scene.
   *
   * @param {plasm.Model}
   * @return {plasm.Scene} for chaining
   * @api public
   */

  plasm.Scene.prototype.add = function (object) {
    if (object instanceof plasm.Camera) {
      this.root.add(object.optics);
      return this;
    }
    if (object instanceof plasm.Model) {
      this.root.add(object.mesh);
      return this;
    }
    return this;
  };

  /**
   * Remove object from the scene.
   *
   * @param {plasm.Model}
   * @return {plasm.Scene} for chaining
   * @api public
   */

  plasm.Scene.prototype.remove = function (object) {
    if (object instanceof plasm.Camera) {
      this.root.remove(object.optics);
      return this;
    }
    if (object instanceof plasm.Model) {
      this.root.remove(object.mesh);
      return this;
    }
    return this;
  };

  /**
   * Return scene bounding sphere's radius.
   *
   * @return {Number} boundingRadius
   * @api public
   */

  plasm.Scene.prototype.getBoundingRadius = function () {
    var radius = 0;
    var geometry;
    var position;
    var maxPos;

    THREE.SceneUtils.traverseHierarchy(this.root, function (obj) {
      if (obj instanceof THREE.Object3D 
          && obj.geometry 
          && obj.geometry.vertices.length) {
        
        geometry = obj.geometry;
        geometry.computeBoundingSphere();
        position = obj.position;
        maxPos = max(abs(position.x), abs(position.y), abs(position.z));
        radius = max(radius, geometry.boundingSphere.radius + maxPos);
      }
    });

    return radius;
  };

  /**
   * Return centroid of the scene.
   *
   * @return {THREE.Vector3} centroid
   * @api public
   */

  plasm.Scene.prototype.getCentroid = function () {
    var centroid = new THREE.Vector3()
    var geometry;
    var position;
    var bbox;
    var p = { 
        max: new THREE.Vector3()
      , min: new THREE.Vector3() 
    };
    
    THREE.SceneUtils.traverseHierarchy(this.root, function (obj) {
      if (obj instanceof THREE.Object3D 
          && obj.geometry 
          && obj.geometry.vertices.length) {

        geometry = obj.geometry;
        geometry.computeBoundingBox();
        bbox = geometry.boundingBox;
        p.max.x = max(p.max.x, bbox.max.x);
        p.max.y = max(p.max.y, bbox.max.y);
        p.max.z = max(p.max.z, bbox.max.z);
        p.min.x = min(p.min.x, bbox.min.x);
        p.min.y = min(p.min.y, bbox.min.y);
        p.min.z = min(p.min.z, bbox.min.z);
      }
    });

    return centroid.add(p.max, p.min).divideScalar(2);
  };

  /**
   * Create a camera.
   *
   * @param {Object} options options
   * @param {Number} [options.fov = 60] field of view
   * @param {Number} [options.aspect = 1] aspect
   * @param {Number} [options.near = 1] near
   * @param {Number} [options.far = 10000] far
   * @return {plasm.Camera} camera
   * @api public
   */

  plasm.Camera = function (options) {
    if (!(this instanceof plasm.Camera)) {
      return new plasm.Camera(options);
    }

    var options = options || {};
    var fovy = options.fovy || 60;
    var aspect = options.aspect || window.innerWidth / window.innerHeight;
    var near = options.near || 0.1;
    var far = options.far || 10000;
    var optics;
    var controls;

    optics = new THREE.PerspectiveCamera(fovy, aspect, near, far);
    optics.fovy = fovy;
    optics.aspect = aspect;
    optics.near = near;
    optics.far = far;
    optics.position.x = 5;
    optics.position.y = 5;
    optics.position.z = 5;

    this.optics = optics;
  };

  /**
   * Create camera controls.
   *
   * @param {plasm.Camera} camera camera to control
   * @param {Object} options options
   * @param {Number} [options.rotateSpeed = 5.0] rotate speed
   * @param {Number} [options.panSpeed = 0.8] pan speed
   * @param {Boolean} [options.noZoom = false] no zoom flag
   * @param {Boolean} [options.noPan = false] no pan flag
   * @param {Boolean} [options.staticMoving = true] static moving flag
   * @param {Number} [options.dynamicDampingFactor = 0.3] dynamic damping factor
   * @param {Array} [options.keys = [65, 83, 68]] control keys
   * @return {plasm.Controls} controls
   * @api public
   */

  plasm.Controls = function (camera, scene, options) {
    if (!(this instanceof plasm.Controls)) {
      return new plasm.Controls(camera, scene, options);
    }

    var options = options || {};
    var controls;

    controls = new THREE.EnhancedTrackballLightControls(camera.optics, scene);
    controls.rotateSpeed = options.rotateSpeed || 5.0;
    controls.zoomSpeed = options.zoomSpeed || 1.2;
    controls.panSpeed = options.panSpeed || 0.8;
    controls.noZoom = options.noZoom || false;
    controls.noPan = options.noPan || false;
    controls.dynamicDampingFactor = options.dynamicDampingFactor || 0.1;

    this.controls = controls;
  };

  /**
   * Update.
   *
   * @return {plasm.Controls} for chaining
   * @api public
   */

  plasm.Controls.prototype.update = function () {
    this.controls.update();
    return this;
  };

  /**
   * Materials
   */

  plasm.materials = {};

  plasm.materials.PointMaterial = function () {
    return new THREE.ParticleBasicMaterial({ 
      color: 0x000000,
      size: 0.2
    });
  };

  plasm.materials.LineMaterial = function () {
    return new THREE.LineBasicMaterial({
        color: 0x292929
      , opacity: 1
      , linewidth: 2
    });
  };

  plasm.materials.MeshMaterial = function () {
    return new THREE.MeshLambertMaterial({
        color: 0xD7D7D7
      , wireframe: !Detector.webgl
      , shading: THREE.FlatShading
      , vertexColors: THREE.FaceColors
    });
  };

  plasm.materials.cloneColor = function (materialFrom, materialTo) {
    var color = materialFrom.color;
    materialTo.color.setRGB(color.r, color.g, color.b);
    materialTo.opacity = materialFrom.opacity !== undefined ? materialFrom.opacity : 1;
    materialTo.transparent = materialFrom.transparent !== undefined ? materialFrom.transparent : false;
  }

  /**
   * Create the model of the given simplicial complex 
   * and add it to the scene of the given viewer.
   *
   * @param {simplexn.SimplicialComplex} complex
   * @param {plasm.Viewer} viewer
   * @return {plasm.Model} model
   * @api public
   */

  exports.DRAW_SINGLE_SIDE = false;

  plasm.Model = function (complex, viewer) {
    if (!(this instanceof plasm.Model)) {
      return new plasm.Model(complex, viewer);
    }

    var complex = complex || new simplexn.SimplicialComplex();
    var pointset = complex.pointset;
    var topology = complex.topology;
    var dim = topology.dim;
    var cells, n_cells, i_cell;
    var v1, v2, v3; 

    var geometry = new THREE.Geometry();
    var material;
    var mesh;

    if (dim <= 0) {
      cells = topology.cells0d();
      n_cells = cells.length;

      for (i_cell = 0; i_cell < n_cells; i_cell += 1) {
        v1 = pointset.get(cells[i_cell]);
        geometry.vertices.push(new THREE.Vertex(
          new THREE.Vector3(v1[0] || 0, v1[1] || 0, v1[2] || 0)
        ));
      }

      material = new plasm.materials.PointMaterial();
      mesh = new THREE.ParticleSystem(geometry, material);
    }

    if (dim === 1) {
      cells = topology.complexes[1];
      n_cells = cells.length;

      for (i_cell = 0; i_cell < n_cells; i_cell += 2) {
        v1 = pointset.get(cells[i_cell + 0]);
        v2 = pointset.get(cells[i_cell + 1]);
        geometry.vertices.push(new THREE.Vertex(
          new THREE.Vector3(v1[0], v1[1], v1[2])
        ));
        geometry.vertices.push(new THREE.Vertex(
          new THREE.Vector3(v2[0], v2[1], v2[2])
        ));
      }

      material = new plasm.materials.LineMaterial();
      mesh = new THREE.Line(geometry, material, THREE.LinePieces);
    }

    if (dim >= 2) {
      cells = topology.complexes[2];
      n_cells = cells.length;

      pointset.forEach(function (v) {
        geometry.vertices.push(new THREE.Vertex(
          new THREE.Vector3(v[0] || 0, v[1] || 0, v[2] || 0)
        ));
      });
    
      for (i_cell = 0; i_cell < n_cells; i_cell += 3) {
        geometry.faces.push(new THREE.Face3(
          cells[i_cell + 0], cells[i_cell + 1], cells[i_cell + 2]
        ));
        if (! exports.DRAW_SINGLE_SIDE) {
          geometry.faces.push(new THREE.Face3(
          cells[i_cell + 2], cells[i_cell + 1], cells[i_cell + 0]
        ));
        }
      }

      geometry.computeCentroids();
      geometry.mergeVertices();
      geometry.computeFaceNormals();

      material = new plasm.materials.MeshMaterial();
      mesh = new THREE.Mesh(geometry, material);
    }

    this.complex = complex;
    this.geometry = geometry;
    this.geometry.dynamic = true;
    this.material = material;
    this.mesh = mesh;
    this.mesh.matrixAutoUpdate = true;
    this.mesh.doubleSided = exports.DRAW_SINGLE_SIDE ? true : false;
    this.viewer = viewer;
  };

  /**
   * Clone
   * 
   * @return {plasm.Clone} clone
   * @api public
   */

  plasm.Model.prototype.clone = function () {
    var model =  new plasm.Model(this.complex.clone(), this.viewer);
    plasm.materials.cloneColor(this.material, model.material);

    return model;
  };

  /**
   * Draw.
   *
   * @return {plasm.Mode} for chaining
   * @api public
   */

  plasm.Model.prototype.draw = function () {
    if (this.mesh.parent !== this.viewer.scene) {
      this.viewer.scene.add(this);
    }

    return this;
  };

  /**
   * Remove this model from the scene.
   *
   * @api public
   */

  plasm.Model.prototype.erase = function () {
    this.viewer.scene.remove(this);
    this.complex = null;
    this.geometry = null;
    this.mesh = null;
  };

  /**
   * Show.
   *
   * @return {plasm.Mode} for chaining
   * @api public
   */

  plasm.Model.prototype.show = function () {
    this.mesh.visible = true;

    return this;
  };

  /**
   * Hide.
   *
   * @return {plasm.Mode} for chaining
   * @api public
   */

  plasm.Model.prototype.hide = function () {
    this.mesh.visible = false;
    
    return this;
  };

  /**
   * Cancel `this` object from graph of the scene.
   *
   * @return {plasm.Model} for chaining
   * @api public
   */

  plasm.Model.prototype.cancel = function () {
    this.mesh.parent.remove(this.mesh);
    return this;
  };

  /**
   * Translate.
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {plasm.Model} this for chaining
   * @api public
   */

  plasm.Model.prototype.translate = function (dims, values) {
    var v = [];
    this.complex.translate(dims, values);

    if (this.complex.rn <= 3) {
      dims.forEach(function (dim, i) {
        v[dim] = values[i];
      });
      this.mesh.position.addSelf({ x: v[0] || 0, y: v[1] || 0, z: v[2] || 0 });
      this.geometry.__dirtyVertices = true;
    }

    return this;
  };

  /**
   * Rotate.
   * 
   * @param {Array|Uint32Array} dims
   * @param {Number|Array|Uint32Array} angle
   * @return {plasm.Model} this for chaining
   * @api public
   */

  plasm.Model.prototype.rotate = function (dims, angle) {
    var v = [];
    var axis = 3 - (dims[0] + dims[1]);
    var angle = angle[0] || angle;

    this.complex.rotate(dims, angle);

    if (this.complex.rn <= 3) {
      v[axis] = angle;
      this.mesh.rotation.addSelf({x: v[0] || 0, y: v[1] || 0, z: v[2] || 0 });
      this.geometry.__dirtyVertices = true;
    }

    return this;
  };

  /**
   * Scale.
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {plasm.Model} this for chaining
   * @api public
   */

  plasm.Model.prototype.scale = function (dims, values) {
    var v = [];
    this.complex.scale(dims, values);

    if (this.complex.rn <= 3) {
      dims.forEach(function (dim, i) {
        v[dim] = values[i];
      });
      this.mesh.scale.multiplySelf({ x: v[0] || 1, y: v[1] || 1, z: v[2] || 1 });
      this.geometry.__dirtyVertices = true;
    }

    return this;
  };

  /**
   * Map.
   * 
   * Map the vertices of this model by the mapping function.
   *
   * @example
   *
   *   domain([[0,1]],[0,2*Math.PI])
   *     .map(function (v) { return [Math.sin(v[0]), Math.cos(v[1])]; })
   *     .draw();
   *
   * @example
   *
   *   domain([[0,1]],[0,2*Math.PI])
   *     .map([
   *        function (v) { return Math.sin(v[0]); }, 
   *        function (v) { return Math.cos(v[1]); }
   *      ])
   *     .draw();
   *
   *
   * @param {Function | Array} mapping
   * @return {plasm.Model} a new mapped model
   * @api public
   */

  plasm.Model.prototype.map = function (mapping, merge) {
    if (mapping instanceof Array) {
      return this.map(function (v) {
        return mapping.map(function (f) {
          return f(v);
        });
      }, merge);
    }

    var complex = this.complex.clone().map(mapping, merge);
    var model = new plasm.Model(complex, this.viewer);

    return model;
  };

  /**
   * Color.
   *
   * @param {Array} rgba rgba
   * @param {Number} [rgba[0] = 0] r
   * @param {Number} [rgba[1] = 0] g
   * @param {Number} [rgba[2] = 0] b
   * @param {Number} [rgba[3] = 1] a
   * @return {plasm.Object} for chaining
   * @api public
   */

  plasm.Model.prototype.color = function (rgba) {
    var a = rgba[3];

    this.material.color.setRGB(rgba[0] || 0, rgba[1] || 0, rgba[2] || 0);
    
    if (a < 1) {
      this.material.opacity = a;
      this.material.transparent = true;
    }
    
    return this;
  };

  /**
   * boundary
   * 
   * @return {plasm.Model} boundary
   * @api public
   */

  plasm.Model.prototype.boundary = function () {
    var complex = this.complex.clone().boundary();
    var boundary = new plasm.Model(complex, this.viewer);
    plasm.materials.cloneColor(this.material, boundary.material);

    return boundary;
  };

  /**
   * skeleton
   * 
   * @param {Number} dim
   * @return {plasm.Model} skeleton
   * @api public
   */

  plasm.Model.prototype.skeleton = function (dim) {
    var complex = this.complex.clone().skeleton(dim);
    var skeleton = new plasm.Model(complex, this.viewer);
    plasm.materials.cloneColor(this.material, skeleton.material);

    return skeleton;
  };

  /**
   * extrude
   * 
   * @param {Array|Float32Array} hlist which must be made by positive numbers 
   *   or by an alternation of positive and negative numbers
   * @return {plasm.Model} extrusion
   * @api private
   */

  plasm.Model.prototype.extrude = function (hlist) {
    var complex = this.complex.clone().extrude(hlist);
    var extrusion = new plasm.Model(complex, this.viewer);
    plasm.materials.cloneColor(this.material, extrusion.material);

    return extrusion;
  };

  /**
   * explode
   *
   * @param {Array|Float32Array} values
   * @return {plasm.Model} explosion
   * @api public
   */

  plasm.Model.prototype.explode = function (values) {
    var complex = this.complex.clone().explode(values);
    var explosion = new plasm.Model(complex, this.viewer);
    plasm.materials.cloneColor(this.material, explosion.material);

    return explosion;
  };

  /**
   * prod
   * At the moment it's customed and tested only for following cases:
   * - 1-rn x 1-rn
   * - 1-rn x 2-rn
   * - 2-rn x 1-rn
   * 
   * @param {plasm.Model} model
   * @return {plasm.Model} result of the product operation
   * @api private
   */

  plasm.Model.prototype._prod = function (model) {
    var complex = this.complex.clone().prod(model.complex);
    var result = new plasm.Model(complex, this.viewer);
    plasm.materials.cloneColor(this.material, result.material);

    return result;
  };

  plasm.Model.prototype.prod1x1 = function (model) {
    if (this.complex.rn !== 1 || model.complex.rn !== 1) throw 'Dimesion error.';

    return this._prod(model);
  };

  plasm.Model.prototype.prod1x2 = function (model) {
    if (this.complex.rn !== 1 || model.complex.rn !== 2) throw 'Dimesion error.';

    return this._prod(model);
  };

  plasm.Model.prototype.prod2x1 = function (model) {
    if (this.complex.rn !== 2 || model.complex.rn !== 1) throw 'Dimesion error.';

    return this._prod(model);
  };

  /**
   * Struct
   *
   * @param {items} complex
   * @param {plasm.Viewer} viewer
   * @return {plasm.Struct} struct
   * @api public
   */

  plasm.Struct = function (items) {
    if (!(this instanceof plasm.Struct)) {
      return new plasm.Struct(items);
    }
    
    var items = items || [];
    var structs = [];
    var models = [];
    var model;

    items.forEach(function (item) {
      model = item.clone();
      if (model instanceof plasm.Model) {
        models.push(model);
      } else if (model instanceof plasm.Struct) {
        structs.push(model);
      }
    });

    this.structs = structs;
    this.models = models;
  };

  /**
   * Clone
   * 
   * @return {plasm.Struct} clone
   * @api public
   */

  plasm.Struct.prototype.clone = function () {
    var cloned = new plasm.Struct();
    var models = [];
    var structs = [];

    this.models.forEach(function (model) {
      models.push(model.clone());
    });

    this.structs.forEach(function (struct) {
      structs.push(struct.clone());
    });

    cloned.models = models;
    cloned.structs = structs;

    return cloned;
  };

  /**
   * Draw.
   *
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.draw = function () {
    this.models.forEach(function (model) {
      model.draw();
    });

    this.structs.forEach(function (struct) {
      struct.draw();
    });

    return this;
  };

  /**
   * Cancel.
   *
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.cancel = function () {
    this.models.forEach(function (model) {
      model.cancel();
    });

    this.structs.forEach(function (struct) {
      struct.cancel();
    });

    return this;
  };

  /**
   * Remove this model from the scene.
   *
   * @api public
   */

  plasm.Struct.prototype.erase = function () {
    this.viewer.scene.remove(this);
    this.structs = null;
    this.models = null;
  };

  /**
   * Show.
   *
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.show = function () {
    this.models.forEach(function (model) {
      model.show();
    });

    this.structs.forEach(function (struct) {
      struct.show();
    });

    return this;
  };

  /**
   * Hide.
   *
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.hide = function () {
    this.models.forEach(function (model) {
      model.hide();
    });

    this.structs.forEach(function (struct) {
      struct.hide();
    });

    return this;
  };

  /**
   * Add `item` to this plasm.Struct.
   *
   * @param {plasm.Model|plasm.Struct} item
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.add = function (item) {
    if (this === item) return this;
    var model = item.clone();
    if (model instanceof plasm.Model) {
      this.models.push(model);
    } else if (model instanceof Struct) {
      this.structs.push(model.complex);
    }

    return this;
  };

  /**
   * Remove `item` from this model.
   *
   * @param {plasm.Model|plasm.Struct} item
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.remove = function (item) {
    if (this === item) return this;
    var index;
    var array;

    if (item instanceof plasm.Model) {
      array = this.models;
    } else if (item instanceof Struct) {
      array = this.structs;
    }
    index = array.indexOf(item);
    if (index > 0) {
      array.splice(index, 1);
      item.remove();
    }
    return this;
  };

  /**
   * Translate.
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {plasm.Struct} this for chaining
   * @api public
   */

  plasm.Struct.prototype.translate = function (dims, values) {
    this.models.forEach(function (model) {
      model.translate(dims, values);
    });

    this.structs.forEach(function (struct) {
      struct.translate(dims, values);
    });

    return this;
  };

  /**
   * Rotate.
   * 
   * @param {Array|Uint32Array} dims
   * @param {Number|Array|Uint32Array} angle
   * @return {plasm.Struct} this for chaining
   * @api public
   */

  plasm.Struct.prototype.rotate = function (dims, angle) {
    this.models.forEach(function (model) {
      model.rotate(dims, angle);
    });

    this.structs.forEach(function (struct) {
      struct.rotate(dims, angle);
    });

    return this;
  };

  /**
   * Scale.
   * 
   * @param {Array|Uint32Array} dims
   * @param {Array|Float32Array} values
   * @return {plasm.Struct} this for chaining
   * @api public
   */

  plasm.Struct.prototype.scale = function (dims, values) {
    this.models.forEach(function (model) {
      model.scale(dims, values);
    });

    this.structs.forEach(function (struct) {
      struct.scale(dims, values);
    });

    return this;
  };

  /**
   * Color.
   *
   * @param {Array} rgba rgba
   * @param {Number} [rgba[0] = 0] r
   * @param {Number} [rgba[1] = 0] g
   * @param {Number} [rgba[2] = 0] b
   * @param {Number} [rgba[3] = 1] a
   * @return {plasm.Struct} for chaining
   * @api public
   */

  plasm.Struct.prototype.color = function (rgba) {
    this.models.forEach(function (model) {
      model.color(rgba);
    });

    this.structs.forEach(function (struct) {
      struct.color(rgba);
    });

    return this;
  };

  /**
   * Extrude.
   *
   * @param {Array|Float32Array} hlist which must be made by positive numbers 
   *   or by an alternation of positive and negative numbers
   * @return {plasm.Struct} extrusion
   * @api private
   */

  plasm.Struct.prototype.extrude = function (hlist) {
    var struct = new plasm.Struct;

    this.models.forEach(function (model) {
      struct.models.push(model.extrude(hlist));
    });

    this.structs.forEach(function (struct) {
      struct.structs.push(struct.extrude(hlist));
    });

    return struct;
  };

  /**
   * SimplicialComplex
   * 
   * @param {Array|Float32Array} points
   * @param {Array|Uint32Array} complexes
   * @return {plasm.Model} simplicial complex
   * @api public
   */
  
  plasm.Viewer.prototype.simplicialComplex = function (points, complex) {
    var complex = new simplexn.SimplicialComplex(points, complex);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * simplex
   * 
   * @param {number} d
   * @return {plasm.Model} a simplex
   * @api public
   */

  plasm.Viewer.prototype.simplex = function (d) {
    var complex = new simplexn.geometries.simplex(d);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * polyline
   * 
   * @param {Array} points
   * @return {simplexn.SimplicialComplex} a polyline
   * @api public
   */

  plasm.Viewer.prototype.polyline = function (points) {
    var complex = simplexn.geometries.polyline(points);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * polypoint
   * 
   * @param {Array} points
   * @return {simplexn.SimplicialComplex} a polypoint
   * @api public
   */

  plasm.Viewer.prototype.polypoint = function (points) {
    var complex = simplexn.geometries.polypoint(points);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * axes
   * 
   * @param {dim} points
   * @faces {Array|Uint32Array} complexes
   * @api public
   */
  
  plasm.Viewer.prototype.axes = function () {
    var axeX = new simplexn.SimplicialComplex([[0,0,0],[1,0,0]],[[0,1]]);
    var axeY = new simplexn.SimplicialComplex([[0,0,0],[0,1,0]],[[0,1]]);
    var axeZ = new simplexn.SimplicialComplex([[0,0,0],[0,0,1]],[[0,1]]);
    var modelX = (new plasm.Model(axeX, this)).color([1,0,0]);
    var modelY = (new plasm.Model(axeY, this)).color([0,1,0]);
    var modelZ = (new plasm.Model(axeZ, this)).color([0,0,1]);
    var axes = new plasm.Struct([modelX,modelY,modelZ]);

    return axes;
  };

  /**
   * simplexGrid
   * 
   * @param {Array} quotesList is a list of hlist made by positive numbers 
   * or made by an alternation of positive and negative numbers
   * @return {plasm.Model} a grid of simplexes
   * @api public
   */

  plasm.Viewer.prototype.simplexGrid = function (quotesList) {
    var complex = new simplexn.geometries.simplexGrid(quotesList);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * cuboid
   * 
   * @param {Array} sides
   * @return {plasm.Model} a cuboidal simplicial complex
   * @api public
   */

  plasm.Viewer.prototype.cuboid = function (sides) {
    var complex = new simplexn.geometries.cuboid(sides);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * intervals
   *
   * @param {Array} values
   * @return {plasm.Model} intervals
   * @api public
   */

  plasm.Viewer.prototype.intervals = function (tip, n) {
    var complex = new simplexn.geometries.intervals(tip, n);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * domain
   *
   * @param {Array} ends
   * @param {Array} ns
   * @return {plasm.Model} domain
   * @api public
   */

  plasm.Viewer.prototype.domain = function (tips, ns) {
    var domain = new simplexn.geometries.domain(tips, ns);
    var model = new plasm.Model(domain, this);
    return model;
  };

  /**
   * cube
   * 
   * @param {Number} dim
   * @return {plasm.Model} a dim-dimendional cube
   * @api public
   */

  plasm.Viewer.prototype.cube = function (d) {
    var complex = new simplexn.geometries.cube(d);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * circle
   * 
   * @param {Number} [radius=1]
   * @param {Number} [n=32] 
   * @return {plasm.Model} a circle
   * @api public
   */

  plasm.Viewer.prototype.circle = function (radius, n) {
    var complex = new simplexn.geometries.circle(radius, n);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * disk
   * 
   * @param {Number} [radius=1]
   * @param {Number} [n=32]
   * @param {Number} [m=1] 
   * @return {plasm.Model} a disk
   * @api public
   */

  plasm.Viewer.prototype.disk = function (radius, n, m) {
    var complex = new simplexn.geometries.disk(radius, n, m);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * cylinderSurface
   * Produces a cylindrical surface of radius r and heigth h.
   * 
   * @param {Number} [r=1]
   * @param {Number} [h=1]
   * @param {Number} [n=16]
   * @param {Number} [m=2] 
   * @return {plasm.Model} a cylindrical surface
   * @api public
   */
   
  plasm.Viewer.prototype.cylinderSurface = function (r, h, n, m) {
    var complex = new simplexn.geometries.cylinderSurface(r, h, n, m);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * torusSurface
   *
   * produces a toroidal surface of radiuses r,R 
   * approximated with n x m x 2 triangles
   *
   * @param {Number} [r_min=1] r_min
   * @param {Number} [r_max=3] r_max
   * @param {Number} [n=12] n
   * @param {Number} [m=8] m
   * @return {plasm.Model} a torus surface
   */

  plasm.Viewer.prototype.torusSurface = function (r_min, r_max, n, m) {
    var complex = new simplexn.geometries.torusSurface(r_min, r_max, n, m);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * torusSolid
   *
   * produces a toroidal surface of radiuses r,R 
   * approximated with n x m x 2 triangles
   *
   * @param {Number} [r=1] r_min
   * @param {Number} [R=3] r_max
   * @param {Number} [n=12] n
   * @param {Number} [m=8] m
   * @param {Number} [p=8] p
   * @return {plasm.Model} a torus solid
   */

  plasm.Viewer.prototype.torusSolid = function (r_min, r_max, n, m, p) {
    var r_min = r_min || 1;
    var r_max = r_max || 3;
    var n = n || 12;
    var m = m || 8;
    var p = p || 8;
    var complex = new simplexn.geometries.torusSolid(r_min, r_max, n, m, p);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * triangleStrip
   * 
   * @param {Array} points
   * @return {plasm.Model} triangle strip
   * @api public
   */

  plasm.Viewer.prototype.triangleStrip = function (points) {
    var complex = new simplexn.geometries.triangleStrip(points);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * triangleFan
   * 
   * @param {Array} points
   * @return {plasm.Model} triangle strip
   * @api public
   */

  plasm.Viewer.prototype.triangleFan = function (points) {
    var complex = new simplexn.geometries.triangleFan(points);
    var model = new plasm.Model(complex, this);
    return model;
  };

  /**
   * helix
   * 
   * @param {Number} [r=1] r
   * @param {Number} [pitch=1] pitch
   * @param {Number} [n=24] n
   * @param {Number} [turns=1] turns
   * @return {plasm.Model} helix
   * @api public   
   */

  plasm.Viewer.prototype.helix = function (r, pitch, n, turns) {
    var complex = new simplexn.geometries.helix(r, pitch, n, turns);
    var model = new plasm.Model(complex, this);
    return model;
  };

}(this));

/*!
 * plasm-fun.js
 * functional environment for plasm.js
 * Copyright (c) 2012 cvd-lab <cvdlab@email.com> (https://github.com/cvd-lab/)
 * MIT License
 */

 !(function (exports) {

  /**
   * Library namespace.
   */

  var fun = exports.fun = {};

  /**
   * Module dependencies.
   */

  Object.keys(ƒ).forEach(function (key) {
    fun[key.toUpperCase()] = ƒ[key];
  });

  /**
   * Library version.
   */

  fun.version = '0.1.5';

  /**
   * Library init.
   */

  var p;

  fun.PLASM = function (viewer) {
    p = viewer;
    fun.globalize();
  };

  /**
   * Library globalization
   */

  fun.globalize = function () {
    fun.globalize = function () {};
    Object.keys(fun).forEach(function (key) {
      exports[key] = fun[key];
    });
  };

  /**
   * DRAW
   * 
   * @param {plasm.Model|plasm.Struct} object
   * @return {plasm.Model|plasm.Struct} object
   * @api public
   */

  fun.DRAW = function (object) {
    if (!(object instanceof plasm.Model) &&
        !(object instanceof plasm.Struct)) {
      return;
    }

    return object.draw();
  };

    /**
   * CANCEL
   * 
   * @param {plasm.Model|plasm.Struct} object
   * @return {plasm.Model|plasm.Struct} object
   * @api public
   */

  fun.CANCEL = function (object) {
    if (!(object instanceof plasm.Model) &&
        !(object instanceof plasm.Struct)) {
      return;
    }

    return object.cancel();
  };

  /**
   * R
   * 
   * @param {Array|Uint32Array} dims
   * @return {Function}
   *   @param {Number|Array|Uint32Array} angle
   *   @return {Function}
   *      @param {plasm.Model|plasm.Struct} object
   *      @return {plasm.Model|plasm.Struct} rotated clone of object
   * @api public
   */

  var R = 
  fun.R = 
  fun.ROTATE = function (dims) {
    return function (angle) {
      return function (object) {
        return object.clone().rotate(dims, angle);
      };
    };
  };

  /**
   * S
   * 
   * @param {Array|Uint32Array} dims
   * @return {Function}
   *   @param {Number} values
   *   @return {Function}
   *      @param {plasm.Model|plasm.Struct} object
   *      @return {plasm.Model|plasm.Struct} scaled clone of object
   * @api public
   */

  var S = 
  fun.S = 
  fun.SCALE = function (dims) {
    return function (values) {
      return function (object) {
        return object.clone().scale(dims, values);
      };
    };
  };

  /**
   * T
   * 
   * @param {Array|Uint32Array} dims
   * @return {Function}
   *   @param {Number} values
   *   @return {Function}
   *      @param {plasm.Model|plasm.Struct} model
   *      @return {plasm.Model|plasm.Struct} translated clone of object
   * @api public
   */

  var T =
  fun.T =
  fun.TRANSLATE = function (dims) {
    return function (values) {
      return function (object) {
       return object.clone().translate(dims, values);
      };
    };
  };

  /**
   * STRUCT
   * 
   * @param {Array} items
   * @return {plasm.Model}
   * @api public
   */

  fun.STRUCT = function (items) {
    var transformations = function (o) {return o;};
    var objects = [];

    temp = [];

    items.forEach(function (item) {
      if (!(item instanceof plasm.Model) && 
          !(item instanceof plasm.Struct)) {
        transformations = COMP2([transformations, item]);
      } else {
        temp.push(APPLY([transformations, item]).clone());
        objects.push(APPLY([transformations, item]));
      }
    });

    return new plasm.Struct(objects, p);
  };

  /**
   * Map.
   * 
   * Map `domain` by `mapping` function.
   *
   * @example
   *
   *   var domain = DOMAIN([[0,1]],[0,2*PI]);
   *   var mapping = function (v) { return [SIN(v[0]), COS(v[1])]; });
   *   var model = MAP(mapping)(domain);
   *   DRAW(model);
   *
   * @example
   *
   *   var domain = DOMAIN([[0,1]],[0,2*PI]);
   *   var mapping = [
   *         function (v) { return SIN(v[0]); }, 
   *         function (v) { return COS(v[1]); }
   *       ]);
   *   var model = MAP(mapping)(domain)
   *   DRAW(model);
   * 
   * @param {Function} mapping
   * @return {Function}
   *   @param {plasm.Model} domain
   *   @return {plasm.Model}
   * @api public
   */

  fun.MAP = function (mapping) {
    return function (domain) {
      return domain.map(mapping);
    };
  };

  /**
   * extrude
   * 
   * @param {Array|Float32Array} hlist a list of positive numbers 
   *   or an alternation of positive and negative numbers
   * @return {Function}
   *   @param {plasm.Model|plasm.Struct} object to extrude
   *   @return {plasm.Model|plasm.Struct} extrusion
   * @api private
   */

  var EXTRUDE =
  fun.EXTRUDE = function (hlist) {
    return function (object) {
      return object.extrude(hlist);
    };
  };

  /**
   * EXPLODE
   *
   * @param {Array|Float32Array} values
   * @return {Function}
   *    @param {plasm.Model} model
   *    @return {plasm.Model} exploded clone of model
   * @api public
   */

  var EXPLODE =
  fun.EXPLODE = function (values) {
    return function (model) {
      return model.explode(values);
    };
  };

  /**
   * SKELETON
   * 
   * @param {Number} dim
   * @return {Function}
   *   @param {plasm.Model} model
   *   @return {plasm.Model} skeleton
   * @api public
   */

  var SKELETON = 
  fun.SKELETON = function (dim) {
    return function (model) {
      return model.skeleton(dim);
    };
  };

  /**
   * BOUNDARY
   * 
   * @param {plasm.Model} model
   * @return {plasm.Model} boundary
   * @api public
   */

  var BOUNDARY = 
  fun.BOUNDARY = function (model) {
    return model.boundary();
  };

  /**
   * Color.
   *
   * @param {Array} rgb rgb
   * @param {Number} [rgb[0] = 0] r
   * @param {Number} [rgb[1] = 0] g
   * @param {Number} [rgb[2] = 0] b
   * @param {Number} [rgb[3] = 1] a
   * @return {Function}
   *   @param {plasm.Model | plasm.Struct} object
   *   @return {plasm.Model | plasm.Struct} colored object
   * @api public
   */

  var COLOR =
  fun.COLOR = function (rgba) {
    return function (object) {
      return object.clone().color(rgba);
    };
  };

  /**
   * SHOW
   *
   * @param {plasm.Model | plasm.Struct} model to show
   * @return {plasm.Model | plasm.Struct} model
   * @api public
   */

  var SHOW = 
  fun.SHOW = function (object) {
    object.show();
  };

  /**
   * HIDE
   *
   * @return {plasm.Model | plasm.Struct} for chaining
   * @api public
   */

  var HIDE = 
  fun.HIDE = function (object) {
    object.hide();
  };

  /**
   * SIMPLICIAL_COMPLEX
   * 
   * @param {Array|Float32Array} points
   * @return {Function}
   *   @param {Array|Uint32Array} cells
   *   @return {plasm.Model} simplicial cells
   * @api public
   */
  
  var SIMPLICIAL_COMPLEX = 
  fun.SIMPLICIAL_COMPLEX = function (points) {
    return function (cells) {
      return p.simplicialComplex(points, cells);
    };
  };

  /**
   * SIMPLEX
   * 
   * @param {number} d
   * @return {plasm.Model} a simplex
   * @api public
   */

  var SIMPLEX =
  fun.SIMPLEX = function (d) {
    return p.simplex(d);
  };

  /**
   * POLYLINE
   * 
   * @param {Array} points
   * @return {plasm.Model} a polyline
   * @api public
   */

  var POLYLINE 
  fun.POLYLINE = function (points) {
    return p.polyline(points);
  };

  /**
   * POLYPOINT
   * 
   * @param {Array} points
   * @return {plasm.Model} a polypoint
   * @api public
   */

  var POLYPOINT
  fun.POLYPOINT = function (points) {
    return p.polypoint(points);
  };

  /**
   * SIMPLEX_GRID
   * 
   * @param {Array} quotesList is a list of hlist made by positive numbers 
   * or made by an alternation of positive and negative numbers
   * @return {plasm.Model} a grid of simplexes
   * @api public
   */

  var SIMPLEX_GRID = 
  fun.SIMPLEX_GRID = function (quotes) {
    return p.simplexGrid(quotes);
  };

  /**
   * CUBE
   * 
   * @param {Number} dim
   * @return {plasm.Model} a dim-dimendional cube
   * @api public
   */

  var CUBE = 
  fun.CUBE = function (d) {
    return p.cube(d);
  };

  /**
   * CUBOID
   * 
   * @param {Array} sides
   * @return {plasm.Model} a cuboidal simplicial complex
   * @api public
   */

  var CUBOID =
  fun.CUBOID = function (sides) {
    return p.cuboid(sides);
  };

  /**
   * INTERVALS
   *
   * @param {Number} tip
   * @return {Function}
   *   @param {Number} n
   *   @return {plasm.Model} intervals
   * @api public
   */

  var INTERVALS = 
  fun.INTERVALS = function (tip) {
    return function (n) {
      return p.intervals(tip, n);
    };
  };

  /**
   * DOMAIN
   *
   * @param {Array} ends
   * @return {Function}
   *   @param {Number} ns
   *   @return {plasm.Model} domain
   * @return {plasm.Model} domain
   * @api public
   */

  var DOMAIN = 
  fun.DOMAIN = function (ends) {
    return function (ns) {
      return p.domain(ends, ns);
    };
  };

  /**
   * PROD
   * cartesian products
   *
   * @param {Array} array
   * @param {plasm.Model} [array[0]] model1
   * @param {plasm.Model} [array[1]] model2
   * @return {plasm.Model} result
   */

  var PROD1x1 = 
  fun.PROD1x1 = function (array) {
    return array[0].prod1x1(array[1]);
  };

  var PROD1x2 = 
  fun.PROD1x2 = function (array) {
    return array[0].prod1x2(array[1]);
  };

  var PROD2x1 = 
  fun.PROD2x1 = function (array) {
    return array[0].prod2x1(array[1]);
  };  

  /**
   * CIRCLE
   * 
   * @param {Number} r radius
   * @return {Function}
   *   @param {Number} n subdivisions
   *   @return {plasm.Model} a circle
   * @api public
   */

  var CIRCLE = 
  fun.CIRCLE = function (r) {
    var r = r || 1;
    return function (n) {
      return p.circle(r, n);
    };
  };

  /**
   * DISK
   * 
   * @param {Number} r radius
   * @return {Function}
   *   @param {Array} divs subdivisions
   *   @param {Number} [divs[0]] slices
   *   @param {Number} [divs[1]] stacks
   *   @return {plasm.Model} a disk
   * @api public
   */

  var DISK = 
  fun.DISK = function (r) {
    var r = r || 1;
    return function (divs) {
      var divs = divs || [];
      var slices = divs[0] || 24;
      var stacks = divs[1] || 3;
      return p.disk(r, slices, stacks);
    };
  };

  /**
   * CYLSURFACE
   * 
   * @param {Array} dims
   * @param {Number} [dims[0]=1] radius
   * @param {Number} [dims[1]=1] height
   * @return {Function}
   *   @param {Array} divs
   *   @param {Number} [divs[0]=16] slices
   *   @param {Number} [divs[1]=2]  stacks
   *   @return {plasm.Model} a cylindrical surface
   * @api public
   */
  
  var CYL_SURFACE = 
  fun.CYL_SURFACE = function (dims) {
    var dims = dims || [];
    var r = dims[0] || 1;
    var h = dims[1] || 1;
    return function (divs) {
      var divs = divs || [];
      var slices = divs[0] || 12;
      var stacks = divs[1] || 8;
      return p.cylinderSurface(r, h, slices, stacks);
    };
  };
 
  /**
   * TORUS_SURFACE
   * 
   * @param {Array} dims
   * @param {Number} [dims[0]=0.1] r min
   * @param {Number} [dims[1]=0.9] r max
   * @return {Function}
   *   @param {Array} divs
   *   @param {Number} [divs[0]=12] slices
   *   @param {Number} [divs[1]=8]  stacks
   *   @return {plasm.Model} a torus surface
   * @api public
   */

  var TORUS_SURFACE = 
  fun.TORUS_SURFACE = function (dims) {
    var dims = dims || [];
    var r_min = dims[0] || 1.0;
    var r_max = dims[1] || 1.9;
    return function (divs) {
      var divs = divs || [];
      var n = divs[0] || 12;
      var m = divs[1] || 12;
      return p.torusSurface(r_min, r_max, n, m);
    };
  };

  /**
   * TORUS_SOLID
   * 
   * @param {Array} dims
   * @param {Number} [dims[0]=0.1] r min
   * @param {Number} [dims[1]=0.9] r max
   * @return {Function}
   *   @param {Array} divs
   *   @param {Number} [divs[0]=12] n
   *   @param {Number} [divs[1]=8]  m
   *   @param {Number} [divs[1]=8]  p
   *   @return {plasm.Model} a torus surface
   * @api public
   */

  var TORUS_SOLID = 
  fun.TORUS_SOLID = function (dims) {
    var dims = dims || [];
    var r_min = dims[0] || 0.1;
    var r_max = dims[1] || 0.9;
    return function (divs) {
      var divs = divs || [];
      var n = divs[0] || 12;
      var m = divs[1] || 8;
      var q = divs[2] || 8;
      return p.torusSolid(r_min, r_max, n, m, q);
    };
  };

  /**
   * TRIANGLE_STRIP
   * 
   * @param {Array} points
   * @return {plasm.Model} triangle strip
   * @api public
   */

  var TRIANGLE_STRIP = 
  fun.TRIANGLE_STRIP = function (points) {
    return p.triangleStrip(points);
  };

  /**
   * TRIANGLEFAN
   * 
   * @param {Array} points
   * @return {plasm.Model} triangle strip
   * @api public
   */

  var TRIANGLE_FAN = 
  fun.TRIANGLE_FAN = function (points) {
    return p.triangleFan(points);
  };

  /**
   * HELIX
   * 
   * @param {Number} [r=1] r
   * @param {Number} [pitch=1] pitch
   * @param {Number} [n=24] n
   * @param {Number} [turns=1] turns
   * @return {plasm.Model} helix
   * @api public   
   */

  var HELIX = 
  fun.HELIX = function (r, pitch, n, turns) {
    return p.helix(r, pitch, n, turns);
  };

  /**
   * CUBIC_HERMITE
   * 
   * @param {Function} sel
   * @return {Function}
   *   @param {Array} args
   *   @return {Function}
   *     @param {Array} point
   *     @return {Funciton}
   * @api public
   */

  var CUBIC_HERMIT = 
  fun.CUBIC_HERMITE = function (sel) {
    return function (args) {
      var p1Fn = args[0];
      var p2Fn = args[1];
      var s1Fn = args[2];
      var s2Fn = args[3];

      return function (point) {
        var u = sel(point);
        var u2 = u * u;
        var u3 = u2 * u;

        var p1 = p1Fn instanceof Function ? p1Fn(point) : p1Fn;
        var p2 = p2Fn instanceof Function ? p2Fn(point) : p2Fn;
        var s1 = s1Fn instanceof Function ? s1Fn(point) : s1Fn;
        var s2 = s2Fn instanceof Function ? s2Fn(point) : s2Fn;

        var rn = p1.length;
        var mapped = new Array(rn);
        var i;

        for (i = 0; i < rn; i += 1) {
          mapped[i] = (2*u3-3*u2+1)*p1[i] + (-2*u3+3*u2)*p2[i]+(u3-2*u2+u)*s1[i]+(u3-u2)*s2[i];
        }

        return mapped;
      };
    };
  };

  /**
   * BEZIER
   * 
   * @param {Function} sel
   * @return {Function}
   *   @param {Array} args
   *   @return {Function}
   *     @param {Array} point
   *     @return {Funciton}
   * @api public
   */

  var BEZIER = 
  fun.BEZIER = function (sel) {  
    return function (args) {
      var n = args.length - 1;
      var controldataFn = args;

      return function (point) {
        var t = sel(point);
        var controldata = new Array(n+1);
        var mapped;
        var rn;
        var weight;
        var crtldata;
        var i, k;

        for (i = 0; i <= n; i += 1) {
          crtldata = controldataFn[i];
          controldata[i] = crtldata instanceof Function ? crtldata(point) : crtldata;
        }

        rn = controldata[0].length;
        mapped = new Array(rn);

        for (i = 0; i < rn; i += 1) {
          mapped[i] = 0.0;
        }

        for (i = 0; i <= n; i += 1) {
          weight = CHOOSE([n,i]) * POW([1-t,n-i]) * POW([t,i]);
          for (k = 0; k < rn; k += 1) {
            mapped[k] += weight * controldata[i][k];
          }
        }
        
        return mapped;
      };
    };
  };

  /**
   * CUBIC_UBSPLINE
   * 
   * @param {Function} domain
   * @return {Function}
   *   @param {Array} args
   * @api public
   */

  var CUBIC_UBSPLINE = 
  fun.CUBIC_UBSPLINE = function (domain) {
    return function (args) {
      var q1Fn = args[0];
      var q2Fn = args[1];
      var q3Fn = args[2];
      var q4Fn = args[3];

      return MAP(function (point) {
        var u = S0(point);
        var u2 = u * u;
        var u3 = u2 * u;
        var rn;
        var mapped;
        var i;

        var q1 = q1Fn instanceof Function ? q1Fn(point) : q1Fn;
        var q2 = q2Fn instanceof Function ? q2Fn(point) : q2Fn;
        var q3 = q3Fn instanceof Function ? q3Fn(point) : q3Fn;
        var q4 = q4Fn instanceof Function ? q4Fn(point) : q4Fn;

        rn = q1.length;
        mapped = new Array(rn);

        for (i = 0; i < rn; i +=1) {
          mapped[i] = (1.0/6.0) * ( (-u3+3*u2-3*u+1)*q1[i] + (3*u3-6*u2+4)*q2[i]+ (-3*u3+3*u2+3*u+1)*q3[i] + (u3)*q4[i]);
        }

        return mapped;
      })(domain);
    };
  };

  /**
   * CUBIC_CARDINAL
   * 
   * @param {Function} domain
   * @param {Number} [h=1]
   * @return {Function}
   *   @param {Array} args
   * @api public
   */

  var CUBIC_CARDINAL =
  fun.CUBIC_CARDINAL = function (domain, h) {
    var h = h !== undefined ? h : 1;

    return function (args) {
      var q1Fn = args[0];
      var q2Fn = args[1];
      var q3Fn = args[2];
      var q4Fn = args[3];

      return MAP(function (point) {
        var u = S0(point);
        var u2 = u * u;
        var u3 = u2 * u;
        var rn;
        var mapped;
        var i;

        var q1 = q1Fn instanceof Function ? q1Fn(point) : q1Fn;
        var q2 = q2Fn instanceof Function ? q2Fn(point) : q2Fn;
        var q3 = q3Fn instanceof Function ? q3Fn(point) : q3Fn;
        var q4 = q4Fn instanceof Function ? q4Fn(point) : q4Fn;

        rn = q1.length;
        mapped = new Array(rn);

        for (i = 0; i < rn; i +=1) {
          mapped[i] = (-h*u3+2*h*u2-h*u)*q1[i] +((2-h)*u3+(h-3)*u2+1)*q2[i] + ((h-2)*u3+(3-2*h)*u2+h*u)*q3[i] + (h*u3-h*u2)*q4[i];
        }

        return mapped;
      })(domain);
    };
  };

  /**
   * SPLINE
   * 
   * @param {Function} curve
   * @return {Function}
   *   @param {Array} points
   *   @return {plasm.Struct}
   * @api public
   */

  var SPLINE =
  fun.SPLINE = function (curve) {
    return function (points) {
      var segments = [];
      var length = points.length;
      var tip = length -4 + 1;
      var slice;
      var i;

      for (i = 0; i < tip; i += 1) {
        slice = points.slice(i,i+4);
        segments.push(curve(slice));
      }
        
      return STRUCT(segments);
    };
  };

  /**
   * DE_BOORD
   * Cox and De Boord coefficients 
   *
   * @api private
   */
  var DE_BOORD = function (T, i, k, t, n) {
    var tmin = T[k-1];
    var tmax = T[n+1];
    var ret, num1, div1, num2, div2;

    // DE_BOORDi1(t)
    if (k === 1) { 
      if ((t >= T[i] && t < T[i+1]) || 
          (t === tmax && t >= T[i] && t <= T[i+1])) {
        return 1;
      } else {
        return 0;
      }
    }

    // DE_BOORDik(t)
    ret = 0;
    num1 = t-T[i];
    div1 = T[i+k-1]-T[i];
    
    if (div1 !== 0) {
      ret += (num1/div1) * DE_BOORD(T,i,k-1,t,n);
    }

    num2 = T[i+k]-t;
    div2 = T[i+k]-T[i+1];
    
    if (div2 !== 0) {
      ret += (num2/div2) * DE_BOORD(T,i+1,k-1,t,n);
    }

    return ret;
  }

  /**
   * BSPLINE
   * 
   * @param {Number} degree
   * @return {Function}
   *   @param {Array} knots
   *   @return {Function}
   *     @param {Array} controls
   *     @return {Funciton}
   * @api public
   */

  var BSPLINE =
  fun.BSPLINE = function (degree) {
    return function (knots) {
      return function (controls) {
        var n = controls.length - 1;
        var m = knots.length -1;
        var k = degree + 1;

        // see http://www.na.iac.cnr.it/~bdv/cagd/spline/B-spline/bspline-curve.html
        if (knots.length !== (n+k+1)) {
          throw "Invalid point/knots/degree for bspline!";
        }

        return function (point) {
          var t = point[0];
          var points = new Array(n);
          var rn;
          var control;
          var mapped;
          var coeff;
          var i, j;

          for (i = 0; i <= n; i += 1) {
            control = controls[i];
            points[i] = control instanceof Function ? control(point) : control;
          }

          rn = points[0].length;
          mapped = new Array(rn);

          for (i = 0; i < rn; i += 1) {
            mapped[i] = 0.0;
          }

          for (i = 0; i <= n; i += 1) {
            coeff = DE_BOORD(knots,i,k,t,n);
            for (j = 0; j < rn; j += 1) {
              mapped[j] += points[i][j] * coeff;
            }
          }

          return mapped;
        };
      };
    };
  };

  /**
   * NUBSPLINE
   * 
   * @param {Number} degree
   * @param {Number} [totpoints=80]
   * @return {Function}
   *   @param {Array} knots
   *   @return {Function}
   *     @param {Array} point
   *     @return {plasm.Model}
   * @api public
   */

  var NUBSPLINE =
  fun.NUBSPLINE = function (degree, totpoints) {
    var totpoints = totpoints !== undefined ? totpoints : 80;
    
    return function (knots) {
      return function (points) {
        var m = knots.length;
        var tmin = SMALLEST(knots);
        var tmax = BIGGEST(knots);
        var tsiz = tmax - tmin;
        var size = totpoints - 1;
        var v = new Array(size + 1);
        var domain;
        var i;

        v[0] = -tmin;
        for (i = 1; i <= size; i += 1) {
          v[i] = tsiz / size;
        }

        domain = SIMPLEX_GRID([v]);
        
        return MAP(BSPLINE(degree)(knots)(points))(domain);
      };
    };
  };

  /**
   * NUBS
   * 
   * @parm {Function} sel
   * @return {Function}
   *   @param {Number} degree
   *   @return {Function}
   *     @param {Array} knots
   *     @return {Function}
   *       @param {Array} controls
   *       @return {Funciton}
   * @api public
   */

  var NUBS =
  fun.NUBS = function (sel) {
    return function (degree) {
      return function (knots) {
        return function (controls) {
          var n = controls.length - 1;
          var knotsLength = knots.length;
          var k = degree + 1;

          // see http://www.na.iac.cnr.it/~bdv/cagd/spline/B-spline/bspline-curve.html
          if (knots.length !== (n+k+1)) {
            throw "Invalid point/knots/degree for bspline!";
          }

          return function (point) {
            var t = sel(point);
            var points = new Array(n);
            var kmin = SMALLEST(knots);
            var kmax = BIGGEST(knots);
            var ksize = kmax - kmin;
            var mappedKnots = new Array(knotsLength);
            var rn;
            var control;
            var mapped;
            var coeff;
            var i, j;

            mappedKnots = knots.map(function (knot) {
              return (knot - kmin) / ksize;
            });

            for (i = 0; i <= n; i += 1) {
              control = controls[i];
              points[i] = control instanceof Function ? control(point) : control;
            }

            rn = points[0].length;
            mapped = new Array(rn);

            for (i = 0; i < rn; i += 1) {
              mapped[i] = 0.0;
            } 

            for (i = 0; i <= n; i += 1) {
              coeff = DE_BOORD(mappedKnots,i,k,t,n);
              for (j = 0; j < rn; j += 1) {
                mapped[j] += points[i][j] * coeff;
              }
            }

            return mapped;
          };
        };
      };
    };
  };

  /**
   * RATIONAL_BSPLINE
   * 
   * @param {Number} degree
   * @return {Function}
   *   @param {Array} knots
   *   @return {Function}
   *     @param {Array} controls
   *     @return {Funciton}
   * @api public
   */

  var RATIONAL_BSPLINE =
  fun.RATIONAL_BSPLINE = function (degree) {
    return function (knots) {
      return function (controls) {
        var bspline = BSPLINE(degree)(knots)(controls);

        return function (point) {
          var mapped = bspline(point);
          var last = mapped.slice(-1)[0];

          // rationalize (== divide for the last value)
          if (last !== 0) {
            mapped = mapped.map(function (item) {
              return item / last;
            });
          }
          
          return mapped.slice(0,-1);
        };
      };
    };
  };

  /**
   * NURBSPLINE
   * 
   * @param {Number} degree
   * @param {Number} [totpoints=80]
   * @return {Function}
   *   @param {Array} knots
   *   @return {Function}
   *     @param {Array} points
   *     @return {plasm.Model}
   * @api public
   */

  var NURBSPLINE =
  fun.NURBSPLINE = function (degree, totpoints) {
    var totpoints = totpoints !== undefined ? totpoints : 80;

    return function (knots) {
      return function (points) {
        var m = knots.length;
        var tmin = SMALLEST(knots);
        var tmax = BIGGEST(knots);
        var tsiz = tmax - tmin;
        var size = totpoints - 1;
        var v = new Array(size + 1);
        var domain;
        var i;

        v[0] = -tmin;
        for (i = 1; i <= size; i += 1) {
          v[i] = tsiz / size;
        }

        domain = SIMPLEX_GRID([v]);
        
        return MAP(RATIONAL_BSPLINE(degree)(knots)(points))(domain);
      };
    };
  };

}(this));