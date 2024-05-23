use crate::{
    samples::{BasePhe, Phe},
    BaseGenot, BaseGenotSnv, GenotSnvRef,
};
//use cmatrix;

/// Missing values are ignroed.
pub fn calc_r2(gsnv_1: &GenotSnvRef, gsnv_2: &GenotSnvRef) -> f64 {
    let count_table_by = count_table_by(gsnv_1, gsnv_2);
    let cnomissing = sum_table3by3(count_table_by);

    let counts_1 = sum_table3by3_row(count_table_by);
    let counts_2 = sum_table3by3_col(count_table_by);

    let mean_1 = mean_3(counts_1, cnomissing);
    let mean_2 = mean_3(counts_2, cnomissing);

    let var_1 = var_3(counts_1, mean_1, cnomissing);
    let var_2 = var_3(counts_2, mean_2, cnomissing);

    let dot_table = ((0, 0, 0), (0, 1, 2), (0, 2, 4));
    let dot = inner_product_3by3(count_table_by, dot_table);

    ((dot as f64) / (cnomissing as f64) - mean_1 * mean_2) / (var_1.sqrt() * var_2.sqrt())
}

fn var_3(counts: (usize, usize, usize), mean: f64, n: usize) -> f64 {
    let square_mean = ((counts.1 + 4 * counts.2) as f64) / (n as f64);
    square_mean - mean * mean
}

fn mean_3(counts: (usize, usize, usize), n: usize) -> f64 {
    ((counts.1 + 2 * counts.2) as f64) / (n as f64)
}

fn sum_tuple3(counts: (usize, usize, usize)) -> usize {
    counts.0 + counts.1 + counts.2
}

fn sum_table3by3_row(counts: Table3By3) -> (usize, usize, usize) {
    (
        sum_tuple3(counts.0),
        sum_tuple3(counts.1),
        sum_tuple3(counts.2),
    )
}

fn sum_table3by3_col(coutns: Table3By3) -> (usize, usize, usize) {
    (
        coutns.0 .0 + coutns.1 .0 + coutns.2 .0,
        coutns.0 .1 + coutns.1 .1 + coutns.2 .1,
        coutns.0 .2 + coutns.1 .2 + coutns.2 .2,
    )
}

fn add_tuple3(sums: (usize, usize, usize), sums_: (usize, usize, usize)) -> (usize, usize, usize) {
    (
        sums.0 + sums_.0,
        sums.1 + sums_.1,
        sums.2 + sums_.2,
        //sums.3 + sums_.3,
    )
}

fn add_table3by3(sums: Table3By3, sums_: Table3By3) -> Table3By3 {
    (
        add_tuple3(sums.0, sums_.0),
        add_tuple3(sums.1, sums_.1),
        add_tuple3(sums.2, sums_.2),
        //add_tuple4(sums.3, sums_.3),
    )
}

// common with cmatrix
// TODO: mv to cmatrix
fn sum_table3by3(counts: Table3By3) -> usize {
    sum_tuple3(counts.0) + sum_tuple3(counts.1) + sum_tuple3(counts.2)
}

// TODO: make struct to imple func -> use Table3By3Ar
type Table3By3 = (
    (usize, usize, usize),
    (usize, usize, usize),
    (usize, usize, usize),
);

fn inner_product_tuple3(counts: (usize, usize, usize), dot_table: (usize, usize, usize)) -> usize {
    counts.0 * dot_table.0 + counts.1 * dot_table.1 + counts.2 * dot_table.2
}

fn inner_product_3by3(counts: Table3By3, dot_table: Table3By3) -> usize {
    inner_product_tuple3(counts.0, dot_table.0)
        + inner_product_tuple3(counts.1, dot_table.1)
        + inner_product_tuple3(counts.2, dot_table.2)
}

// ignore missing
fn add_count_3by3_interaction(table: Table3By3, (g1, g2): (u8, u8)) -> Table3By3 {
    let mut table = table;
    match (g1, g2) {
        (0, 0) => table.0 .0 = table.0 .0 + 1,
        (0, 1) => table.0 .1 = table.0 .1 + 1,
        (0, 2) => table.0 .2 = table.0 .2 + 1,
        (1, 0) => table.1 .0 = table.1 .0 + 1,
        (1, 1) => table.1 .1 = table.1 .1 + 1,
        (1, 2) => table.1 .2 = table.1 .2 + 1,
        (2, 0) => table.2 .0 = table.2 .0 + 1,
        (2, 1) => table.2 .1 = table.2 .1 + 1,
        (2, 2) => table.2 .2 = table.2 .2 + 1,
        _ => (), // missing
                 // _ => panic!("Invalid genotype."),
    };
    table
}

//fn calc_mean(gsnv: &GenotSnvRef) -> f64 {
//    gsnv.maf() * 2.0
//}

//fn calc_mean_var(gsnv: &GenotSnvRef) -> f64 {
//    let (c0, c1, c2, _cm) = gsnv.count_table();
//
//    //if cm!=0{
//    //    panic!("Missing values are not allowed.");
//    //}
//
//    //let n=gnsv.n();
//
//    let cnomissing = c0 + c1 + c2;
//    //let cnomissing = n - cm;
//
//    let mean = ((c1 + 2 * c2) as f64) / (cnomissing as f64);
//    let square_mean = ((c1 + 4 * c2) as f64) / (cnomissing as f64);
//
//    let var = square_mean - mean * mean;
//    //let var = square_mean - mean;
//
//    var
//}

type Table4By4 = (
    (usize, usize, usize, usize),
    (usize, usize, usize, usize),
    (usize, usize, usize, usize),
    (usize, usize, usize, usize),
);

fn add_tuple4(
    sums: (usize, usize, usize, usize),
    sums_: (usize, usize, usize, usize),
) -> (usize, usize, usize, usize) {
    (
        sums.0 + sums_.0,
        sums.1 + sums_.1,
        sums.2 + sums_.2,
        sums.3 + sums_.3,
    )
}

#[allow(dead_code)]
fn add_table4by4(sums: Table4By4, sums_: Table4By4) -> Table4By4 {
    (
        add_tuple4(sums.0, sums_.0),
        add_tuple4(sums.1, sums_.1),
        add_tuple4(sums.2, sums_.2),
        add_tuple4(sums.3, sums_.3),
    )
}

#[allow(dead_code)]
fn sum_tuple4(sums: (usize, usize, usize, usize)) -> usize {
    sums.0 + sums.1 + sums.2 + sums.3
}

#[allow(dead_code)]
fn sum_table4by4(sums: Table4By4) -> usize {
    sum_tuple4(sums.0) + sum_tuple4(sums.1) + sum_tuple4(sums.2) + sum_tuple4(sums.3)
}

// TODO: test
// TODO: use u32 for speed ???
// For g1, g2 interaction
// ((x00, x01, x02),(x10,x11,x12),(x20,x21,x22)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Table3By3Ar([usize; 9]);

impl Table3By3Ar {
    pub fn inner(self) -> [usize; 9] {
        self.0
    }

    pub fn new_empty() -> Self {
        Table3By3Ar([0usize; 9])
    }

    pub fn val(&self, g1: usize, g2: usize) -> usize {
        if g1 >= 3 || g2 >= 3 {
            panic!("Invalid genotype.");
        }
        self.0[g1 * 3 + g2]
    }

    pub fn to_3x3(self) -> [[usize; 3]; 3] {
        let x = self.0;
        //let Table3by3Ar(x) = self;
        [[x[0], x[1], x[2]], [x[3], x[4], x[5]], [x[6], x[7], x[8]]]
    }

    pub fn sum(self) -> usize {
        self.0.iter().sum()
    }

    pub fn add(self, other: Self) -> Self {
        let mut x = self.0;
        let y = other.0;
        for i in 0..9 {
            x[i] += y[i];
        }
        Self(x)
    }

    pub fn subtract(self, other: Self) -> Self {
        let mut res = [0usize; 9];
        for i in 0..9 {
            res[i] = self.0[i] - other.0[i];
        }
        Self(res)
    }

    pub fn add_count(self, (g1, g2): (u8, u8)) -> Table3By3Ar {
        let mut x = self.0;
        if g1 >= 3 || g2 >= 3 {
            // do nth
            //panic!("Invalid genotype.");
            return self;
        }
        // 0<=g1, g2<=2
        x[(g1 * 3 + g2) as usize] += 1;
        Table3By3Ar(x)
    }
}

// do not use simd for ni=n/32
// should be faster than count_table_by_slow
fn count_table_by(gsnv_1: &GenotSnvRef, gsnv_2: &GenotSnvRef) -> Table3By3 {
    let n = gsnv_1.n();
    let pred_1_s0m = gsnv_1.predict_s(0);
    let pred_1_s1m = gsnv_1.predict_s(1);
    let pred_2_s0m = gsnv_2.predict_s(0);
    let pred_2_s1m = gsnv_2.predict_s(1);

    let mut sums = ((0, 0, 0), (0, 0, 0), (0, 0, 0));
    //let mut sums = ((0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0));

    // assume no padding here
    //fn sum_byte_count(p1_0: u32, p1_1: u32, p2_0: u32, p2_1: u32) -> (usize, usize, usize) {
    fn sum_byte_count(p1_0: u32, p1_1: u32, p2_0: u32, p2_1: u32) -> Table3By3 {
        // do not use (0,0) since could be mixed with padding
        let v1_2 = p1_0 & p1_1;
        let v1_1 = p1_0 & (!p1_1);
        let v1_0 = (!p1_0) & (!p1_1);
        //let v1_m = (!p1_0) & p1_1;
        let v2_2 = p2_0 & p2_1;
        let v2_1 = p2_0 & (!p2_1);
        let v2_0 = (!p2_0) & (!p2_1);
        //let v2_m = (!p2_0) & p2_1;

        (
            (
                cmatrix::popcnt(v1_0 & v2_0),
                cmatrix::popcnt(v1_0 & v2_1),
                cmatrix::popcnt(v1_0 & v2_2),
                //cmatrix::popcnt(v1_0 & v2_m),
            ),
            (
                cmatrix::popcnt(v1_1 & v2_0),
                cmatrix::popcnt(v1_1 & v2_1),
                cmatrix::popcnt(v1_1 & v2_2),
                //cmatrix::popcnt(v1_1 & v2_m),
            ),
            (
                cmatrix::popcnt(v1_2 & v2_0),
                cmatrix::popcnt(v1_2 & v2_1),
                cmatrix::popcnt(v1_2 & v2_2),
                //cmatrix::popcnt(v1_2 & v2_m),
            ),
            //(
            //    cmatrix::popcnt(v1_m & v2_0),
            //    cmatrix::popcnt(v1_m & v2_1),
            //    cmatrix::popcnt(v1_m & v2_2),
            //    cmatrix::popcnt(v1_m & v2_m),
            //),
        )
    }

    //for ni in 0..(n / 32 + 1) {
    for ni in 0..(n / 32) {
        let pred_1_s0_b32 =
            u32::from_le_bytes(pred_1_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_1_s1_b32 =
            u32::from_le_bytes(pred_1_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let pred_2_s0_b32 =
            u32::from_le_bytes(pred_2_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_2_s1_b32 =
            u32::from_le_bytes(pred_2_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let sums_32 = sum_byte_count(pred_1_s0_b32, pred_1_s1_b32, pred_2_s0_b32, pred_2_s1_b32);

        //sums = add_tuple3(sums, sums_32);
        sums = add_table3by3(sums, sums_32);
    }
    // for ni=n/32
    for sample_i in (n / 32) * 32..n {
        //let pred_1_s0_b32 = gsnv_1.predict_s(0)[sample_i];
        // let pred_1_s1_b32 = gsnv_1.predict_s(1)[sample_i];
        // let pred_2_s0_b32 = gsnv_2.predict_s(0)[sample_i];
        // let pred_2_s1_b32 = gsnv_2.predict_s(1)[sample_i];

        let g1 = gsnv_1.get_val_unchecked(sample_i);
        let g2 = gsnv_2.get_val_unchecked(sample_i);

        //let sums_32 = sum_byte_count(pred_1_s0_b32, pred_1_s1_b32, pred_2_s0_b32, pred_2_s1_b32);

        sums = add_count_3by3_interaction(sums, (g1, g2));
    }

    //println!("in table afr for: {} sec",  start_time.elapsed().as_micros());

    //let c2 = sums.0;
    //let c1 = sums.1;
    //let cm = sums.2;

    // let c_exclude_0by0 = sum_table4by4(sums);
    // let c_0by0 = n - c_exclude_0by0;

    // sums.0 .0 = c_0by0;

    //let c0 = call - (c2 + c1 + cm);

    //(c0, c1, c2, cm)

    // let sums_3by3 = (
    //     (sums.0 .0, sums.0 .1, sums.0 .2),
    //     (sums.1 .0, sums.1 .1, sums.1 .2),
    //     (sums.2 .0, sums.2 .1, sums.2 .2),
    // );

    sums
}

// for test
#[allow(dead_code)]
fn count_table_by_phe_naive(
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
    phe: &Phe,
) -> (Table3By3Ar, Table3By3Ar, Table3By3Ar) {
    let n = gsnv_1.n();
    assert_eq!(n, gsnv_2.n());

    let mut sums_case = Table3By3Ar::new_empty();
    let mut sums_cont = Table3By3Ar::new_empty();

    for sample_i in 0..n {
        let g1 = gsnv_1.get_val_unchecked(sample_i);
        let g2 = gsnv_2.get_val_unchecked(sample_i);
        let y = phe.get_unchecked(sample_i);

        if y {
            sums_case = sums_case.add_count((g1, g2));
        } else {
            sums_cont = sums_cont.add_count((g1, g2));
        }
    }

    let sums_both = sums_case.add(sums_cont);

    (sums_case, sums_cont, sums_both)
}

// do not use simd for ni=n/32
// should be faster than count_table_by_slow
// return table for (case, control, both)
pub fn count_table_by_phe(
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
    phe: &Phe,
) -> (Table3By3Ar, Table3By3Ar, Table3By3Ar) {
    let n = gsnv_1.n();
    let pred_1_s0m = gsnv_1.predict_s(0);
    let pred_1_s1m = gsnv_1.predict_s(1);
    let pred_2_s0m = gsnv_2.predict_s(0);
    let pred_2_s1m = gsnv_2.predict_s(1);
    let ys = phe.inner();

    let mut sums_both = Table3By3Ar::new_empty();
    let mut sums_case = Table3By3Ar::new_empty();
    //let mut sums_both = ((0, 0, 0), (0, 0, 0), (0, 0, 0));
    //let mut sums_case = ((0, 0, 0), (0, 0, 0), (0, 0, 0));

    // assume no padding here
    fn sum_byte_count(
        y: u32,
        p1_0: u32,
        p1_1: u32,
        p2_0: u32,
        p2_1: u32,
    ) -> (Table3By3Ar, Table3By3Ar) {
        // you can use (0,0) since no padding until ni<n/32
        let v1_2 = p1_0 & p1_1;
        let v1_1 = p1_0 & (!p1_1);
        let v1_0 = (!p1_0) & (!p1_1);
        let v2_2 = p2_0 & p2_1;
        let v2_1 = p2_0 & (!p2_1);
        let v2_0 = (!p2_0) & (!p2_1);

        let vi_00 = v1_0 & v2_0;
        let vi_01 = v1_0 & v2_1;
        let vi_02 = v1_0 & v2_2;
        let vi_10 = v1_1 & v2_0;
        let vi_11 = v1_1 & v2_1;
        let vi_12 = v1_1 & v2_2;
        let vi_20 = v1_2 & v2_0;
        let vi_21 = v1_2 & v2_1;
        let vi_22 = v1_2 & v2_2;

        (
            Table3By3Ar([
                cmatrix::popcnt(vi_00),
                cmatrix::popcnt(vi_01),
                cmatrix::popcnt(vi_02),
                cmatrix::popcnt(vi_10),
                cmatrix::popcnt(vi_11),
                cmatrix::popcnt(vi_12),
                cmatrix::popcnt(vi_20),
                cmatrix::popcnt(vi_21),
                cmatrix::popcnt(vi_22),
            ]),
            Table3By3Ar([
                cmatrix::popcnt(y & vi_00),
                cmatrix::popcnt(y & vi_01),
                cmatrix::popcnt(y & vi_02),
                cmatrix::popcnt(y & vi_10),
                cmatrix::popcnt(y & vi_11),
                cmatrix::popcnt(y & vi_12),
                cmatrix::popcnt(y & vi_20),
                cmatrix::popcnt(y & vi_21),
                cmatrix::popcnt(y & vi_22),
            ]),
        )
    }

    //for ni in 0..(n / 32 + 1) {
    for ni in 0..(n / 32) {
        let pred_1_s0_b32 =
            u32::from_le_bytes(pred_1_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_1_s1_b32 =
            u32::from_le_bytes(pred_1_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let pred_2_s0_b32 =
            u32::from_le_bytes(pred_2_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_2_s1_b32 =
            u32::from_le_bytes(pred_2_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let (sums_32_both, sums_32_case) = sum_byte_count(
            ys_b32,
            pred_1_s0_b32,
            pred_1_s1_b32,
            pred_2_s0_b32,
            pred_2_s1_b32,
        );

        sums_both = sums_both.add(sums_32_both);
        sums_case = sums_case.add(sums_32_case);
        //sums_both = add_table3by3(sums_both, sums_32_both);
        //sums_case = add_table3by3(sums_case, sums_32_case);
        //sums = add_tuple3(sums, sums_32);
    }
    // for ni=n/32
    for sample_i in (n / 32) * 32..n {
        let g1 = gsnv_1.get_val_unchecked(sample_i);
        let g2 = gsnv_2.get_val_unchecked(sample_i);
        let y = phe.get_unchecked(sample_i);

        sums_both = sums_both.add_count((g1, g2));
        //sums_both = add_count_3by3_interaction(sums_both, (g1, g2));
        if y {
            sums_case = sums_case.add_count((g1, g2));
            //sums_case = add_count_3by3_interaction(sums_case, (g1, g2));
        }
    }

    let sums_cont = sums_both.subtract(sums_case);

    (sums_case, sums_cont, sums_both)
    //(sums_both,sums_case, sums_cont)
    //(sums_case, sums_cont)
}

// for test
#[allow(dead_code)]
fn count_table_by_old(gsnv_1: &GenotSnvRef, gsnv_2: &GenotSnvRef) -> Table3By3 {
    let n = gsnv_1.n();
    let pred_1_s0m = gsnv_1.predict_s(0);
    let pred_1_s1m = gsnv_1.predict_s(1);
    let pred_2_s0m = gsnv_2.predict_s(0);
    let pred_2_s1m = gsnv_2.predict_s(1);

    //let mut sums = ((0, 0, 0), (0, 0, 0), (0, 0, 0));
    let mut sums = ((0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0));

    //fn sum_byte_count(p1_0: u32, p1_1: u32, p2_0: u32, p2_1: u32) -> (usize, usize, usize) {
    fn sum_byte_count(p1_0: u32, p1_1: u32, p2_0: u32, p2_1: u32) -> Table4By4 {
        // do not use (0,0) since could be mixed with padding
        let v1_2 = p1_0 & p1_1;
        let v1_1 = p1_0 & (!p1_1);
        let v1_0 = (!p1_0) & (!p1_1);
        let v1_m = (!p1_0) & p1_1;
        let v2_2 = p2_0 & p2_1;
        let v2_1 = p2_0 & (!p2_1);
        let v2_0 = (!p2_0) & (!p2_1);
        let v2_m = (!p2_0) & p2_1;

        (
            (
                0,
                cmatrix::popcnt(v1_0 & v2_1),
                cmatrix::popcnt(v1_0 & v2_2),
                cmatrix::popcnt(v1_0 & v2_m),
            ),
            (
                cmatrix::popcnt(v1_1 & v2_0),
                cmatrix::popcnt(v1_1 & v2_1),
                cmatrix::popcnt(v1_1 & v2_2),
                cmatrix::popcnt(v1_1 & v2_m),
            ),
            (
                cmatrix::popcnt(v1_2 & v2_0),
                cmatrix::popcnt(v1_2 & v2_1),
                cmatrix::popcnt(v1_2 & v2_2),
                cmatrix::popcnt(v1_2 & v2_m),
            ),
            (
                cmatrix::popcnt(v1_m & v2_0),
                cmatrix::popcnt(v1_m & v2_1),
                cmatrix::popcnt(v1_m & v2_2),
                cmatrix::popcnt(v1_m & v2_m),
            ),
        )

        //(
        //    cmatrix::popcnt(c2),
        //    cmatrix::popcnt(c1),
        //    cmatrix::popcnt(cm),
        //)
    }

    //fn add_tuple3(
    //    sums: (usize, usize, usize),
    //    sums_: (usize, usize, usize),
    //) -> (usize, usize, usize) {
    //    (sums.0 + sums_.0, sums.1 + sums_.1, sums.2 + sums_.2)
    //}

    for ni in 0..(n / 32 + 1) {
        let pred_1_s0_b32 =
            u32::from_le_bytes(pred_1_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_1_s1_b32 =
            u32::from_le_bytes(pred_1_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let pred_2_s0_b32 =
            u32::from_le_bytes(pred_2_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_2_s1_b32 =
            u32::from_le_bytes(pred_2_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let sums_32 = sum_byte_count(pred_1_s0_b32, pred_1_s1_b32, pred_2_s0_b32, pred_2_s1_b32);

        //sums = add_tuple3(sums, sums_32);
        sums = add_table4by4(sums, sums_32);
    }

    //println!("in table afr for: {} sec",  start_time.elapsed().as_micros());

    //let c2 = sums.0;
    //let c1 = sums.1;
    //let cm = sums.2;

    let c_exclude_0by0 = sum_table4by4(sums);
    let c_0by0 = n - c_exclude_0by0;

    sums.0 .0 = c_0by0;

    //let c0 = call - (c2 + c1 + cm);

    //(c0, c1, c2, cm)

    let sums_3by3 = (
        (sums.0 .0, sums.0 .1, sums.0 .2),
        (sums.1 .0, sums.1 .1, sums.1 .2),
        (sums.2 .0, sums.2 .1, sums.2 .2),
    );

    sums_3by3
}

//// old
//fn calculate_r2_nosimd(
//    predicts: &[B8],
//    predicts_chosen: &[B8],
//    gt_mean: f64,
//    gt_var: f64,
//    gt_mean_chosen: f64,
//    gt_var_chosen: f64,
//    n: usize,
//) -> f64 {
//    let len_n = genot_twin::len_n(n);
//    //let len_n = genotype::len_n(n);
//    let preds0 = &predicts[..len_n];
//    let preds1 = &predicts[len_n..];
//    let preds0_c = &predicts_chosen[..len_n];
//    let preds1_c = &predicts_chosen[len_n..];
//
//    //let mut dot = 0.0f64;
//    let mut dot: usize = 0;
//    for ni in 0..n {
//        let g = convert_predicts_to_genotype(bit::bget(preds0, ni), bit::bget(preds1, ni));
//        let g_c = convert_predicts_to_genotype(bit::bget(preds0_c, ni), bit::bget(preds1_c, ni));
//        dot += (g as usize) * (g_c as usize);
//        //dot += (g as f64) * (g_c as f64);
//    }
//
//    let r2 = ((dot as f64) / (n as f64) - gt_mean * gt_mean_chosen)
//        / (gt_var.sqrt() * gt_var_chosen.sqrt());
//    r2
//}

#[cfg(test)]
mod tests {
    use super::*;

    //use super::mistake;
    //use super::genotype;
    //use genetics::{plink, SnvId};

    //fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
    //    (v - w).abs() < e
    //}

    //fn setup_test() -> (Vec<B8>, usize, usize) {
    //    let fin = String::from("../../test/data/toy1/genot");
    //    let fin_snv = None;
    //    let fin_sample = None;

    //    let m_in: usize = plink::compute_num_snv(&fin).unwrap();
    //    println!("m_in: {}", m_in);
    //    let n_in: usize = plink::compute_num_sample(&fin).unwrap();
    //    println!("n_in: {}", n_in);
    //    // load snvs
    //    let snvs_in: Vec<Snv> = plink::load_snvs(&fin, m_in);
    //    let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
    //    //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
    //    let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
    //    //let ys: Vec<B8> = plink::load_ys_b8(&fin, n, &use_samples);

    //    //(fin, ys, m, n, use_snvs, use_samples)

    //    let predictions = predict::generate_predictions(&fin, m, n, &use_snvs, &use_samples);

    //    (predictions, m, n)
    //}

    use crate::GenotSnv;

    fn setup_test_2() -> (GenotSnv, GenotSnv) {
        //let (m, n) = (2, 4);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];

        let x0: [u8; 4] = [0, 1, 1, 2];
        let x1: [u8; 4] = [2, 1, 1, 0];

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        //let mistakes = mistake::generate_mistakes_from_x_snv(&x, &ys);
        //let predict0 = predict::generate_predictions_from_x_snv(&x0);
        //let predict1 = predict::generate_predictions_from_x_snv(&x1);
        //println!("predict0: {:?}", predict0);

        //(predict0, predict1, m, n)
        (gsnv0, gsnv1)
    }

    fn setup_test_3() -> (GenotSnv, GenotSnv) {
        // larger test for n>8
        //let (m, n) = (2, 9);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let x0: [u8; 9] = [0, 1, 1, 2, 2, 2, 1, 0, 2];
        let x1: [u8; 9] = [2, 0, 1, 0, 1, 2, 0, 2, 1];

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        (gsnv0, gsnv1)
    }

    // test corner len=33 for simd
    fn setup_test_4() -> (GenotSnv, GenotSnv) {
        // larger test for n>8
        //let (m, n) = (2, 32);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let mut x0 = Vec::with_capacity(33);
        // cannot be all the same values
        x0.push(1u8);
        for _ in 0..32 {
            x0.push(0u8);
        }
        let x1 = x0.clone();

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        (gsnv0, gsnv1)
    }

    // when all vals are 0
    // r2 should be nan
    //fn setup_test_5() -> (Vec<B8>, Vec<B8>, usize, usize) {
    fn setup_test_5() -> (GenotSnv, GenotSnv) {
        //let (m, n) = (2, 4);
        let x0: [u8; 4] = [0, 0, 0, 0];
        let x1: [u8; 4] = [2, 1, 1, 0];

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        (gsnv0, gsnv1)
    }

    // when all vals are the same 1
    // r2 should be nan
    fn setup_test_6() -> (GenotSnv, GenotSnv) {
        let x0: [u8; 4] = [1, 1, 1, 1];
        let x1: [u8; 4] = [2, 1, 1, 0];

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        (gsnv0, gsnv1)
    }

    // includes (0,0)
    fn setup_test_7() -> (GenotSnv, GenotSnv) {
        // larger test for n>8
        //let (m, n) = (2, 9);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let x0 = vec![0u8, 1, 1, 2, 0, 2, 1, 0, 2, 3, 2];
        let x1 = vec![0u8, 0, 1, 0, 0, 2, 0, 2, 1, 0, 3];

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        (gsnv0, gsnv1)
    }

    fn setup_test_7_missing() -> (GenotSnv, GenotSnv) {
        // larger test for n>8
        //let (m, n) = (2, 9);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let x0 = vec![0u8, 1, 1, 2, 0, 2, 1, 0, 2];
        let x1 = vec![0u8, 0, 1, 0, 0, 2, 0, 2, 1];

        let gsnv0 = GenotSnv::new(&x0);
        let gsnv1 = GenotSnv::new(&x1);

        (gsnv0, gsnv1)
    }

    fn setup_test_8() -> (GenotSnv, GenotSnv, Phe) {
        // large test
        let mut vec1 = vec![0u8; 30];
        let vec1_1 = vec![1u8; 30];
        let vec1_2 = vec![2u8; 30];
        vec1.extend(vec1_1);
        vec1.extend(vec1_2);

        let vec2 = vec![vec![0, 1, 2]; 30];
        let vec2 = vec2.into_iter().flatten().collect::<Vec<u8>>();

        let gsnv0 = GenotSnv::new(&vec1);
        let gsnv1 = GenotSnv::new(&vec2);

        let ys = vec![vec![false, true]; 45];
        let ys = ys.into_iter().flatten().collect::<Vec<bool>>();
        let phe = Phe::new(&ys);

        (gsnv0, gsnv1, phe)
    }

    fn setup_test_9_nopad() -> (GenotSnv, GenotSnv, Phe) {
        // large test len=64
        let mut vec1 = vec![0u8; 22];
        let vec1_1 = vec![1u8; 22];
        let vec1_2 = vec![2u8; 20];
        vec1.extend(vec1_1);
        vec1.extend(vec1_2);

        let vec2 = vec![vec![0, 1, 2]; 22];
        let vec2 = vec2.into_iter().flatten().collect::<Vec<u8>>();
        let vec2 = vec2[..64].to_vec();

        let gsnv0 = GenotSnv::new(&vec1);
        let gsnv1 = GenotSnv::new(&vec2);

        let ys = vec![vec![false, true]; 32];
        let ys = ys.into_iter().flatten().collect::<Vec<bool>>();
        let phe = Phe::new(&ys);

        (gsnv0, gsnv1, phe)
    }

    //#[test]
    //fn test_calculate_mean_m_2() {
    //    //let x0: [u8; 4] = [0, 1, 1, 2]; mean: 1, var: 0.5
    //    //let x1: [u8; 4] = [2, 1, 1, 0]; mean: 1, var: 0.5
    //    let mean0_exp = 1.0;
    //    let mean1_exp = 1.0;

    //    let (predict0, predict1) = setup_test_2();

    //    let mean0 = calculate_mean_m(&predict0, n);
    //    let mean1 = calculate_mean_m(&predict1, n);

    //    println!("mean0 {}", mean0);
    //    assert!(is_eq_f64(mean0, mean0_exp, 1e-7));
    //    assert!(is_eq_f64(mean1, mean1_exp, 1e-7));
    //}
    //#[test]
    //fn test_calculate_var_m_2() {
    //    //let x0: [u8; 4] = [0, 1, 1, 2]; mean: 1, var: 0.5
    //    //let x1: [u8; 4] = [2, 1, 1, 0]; mean: 1, var: 0.5
    //    let var0_exp = 0.5;
    //    let var1_exp = 0.5;

    //    let (predict0, predict1, _, n) = setup_test_2();

    //    let mean0 = calculate_mean_m(&predict0, n);
    //    let var0 = calculate_var_m(mean0, &predict0, n);
    //    let mean1 = calculate_mean_m(&predict1, n);
    //    let var1 = calculate_var_m(mean1, &predict1, n);

    //    println!("var0 {}", var0);
    //    assert!(is_eq_f64(var0, var0_exp, 1e-7));
    //    assert!(is_eq_f64(var1, var1_exp, 1e-7));
    //}

    //#[test]
    //fn test_calculate_r2_simd() {
    //    let (predictions, _, n) = setup_test();

    //    let predicts = predict::predictions_snv_s(&predictions, 0, n);
    //    let predicts_chosen = predict::predictions_snv_s(&predictions, 1, n);

    //    let mean = calculate_mean_m(predicts, n);
    //    let var = calculate_var_m(mean, predicts, n);
    //    let mean_c = calculate_mean_m(predicts_chosen, n);
    //    let var_c = calculate_var_m(mean_c, predicts_chosen, n);

    //    unsafe {
    //        let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

    //        let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
    //        println!("r2 {}", r2);

    //        let r2_nosimd =
    //            calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);

    //        assert!(is_eq_f64(r2, r2_nosimd, 1e-7));
    //    }
    //}

    //#[test]
    //fn test_calculate_r2_nosimd() {
    //    let (predictions, _, n) = setup_test();

    //    let predicts = predict::predictions_snv_s(&predictions, 0, n);
    //    let predicts_chosen = predict::predictions_snv_s(&predictions, 1, n);

    //    let mean = calculate_mean_m(predicts, n);
    //    let var = calculate_var_m(mean, predicts, n);
    //    let mean_c = calculate_mean_m(predicts_chosen, n);
    //    let var_c = calculate_var_m(mean_c, predicts_chosen, n);

    //    let r2 = calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);
    //    println!("r2 {}", r2);
    //}

    //#[test]
    //fn test_calculate_r2_self_simd() {
    //    let (predictions, _, n) = setup_test();

    //    let predicts = predict::predictions_snv_s(&predictions, 0, n);
    //    let predicts_chosen = predict::predictions_snv_s(&predictions, 0, n);

    //    let mean = calculate_mean_m(predicts, n);
    //    let var = calculate_var_m(mean, predicts, n);
    //    let mean_c = calculate_mean_m(predicts_chosen, n);
    //    let var_c = calculate_var_m(mean_c, predicts_chosen, n);

    //    unsafe {
    //        let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

    //        let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
    //        println!("r2 {}", r2);

    //        assert!(is_eq_f64(r2, 1.0, 1e-7));
    //    }
    //}

    #[test]
    fn test_count_table_by() {
        let (gsnv1, gsnv2) = setup_test_3();
        let count_table = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        let count_table_exp = ((0, 0, 2), (2, 1, 0), (1, 2, 1));
        assert_eq!(count_table, count_table_exp);
    }

    #[test]
    fn test_count_table_by_7() {
        let (gsnv1, gsnv2) = setup_test_7();
        let count_table = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        let count_table_exp = ((2, 0, 1), (2, 1, 0), (1, 1, 1));
        assert_eq!(count_table, count_table_exp);

        // missing
        let (gsnv1, gsnv2) = setup_test_7_missing();
        let count_table = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        assert_eq!(count_table, count_table_exp);
    }

    #[test]
    fn test_count_table_by_vs_old_7() {
        let (gsnv1, gsnv2) = setup_test_7();
        let count_table = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        let count_table_old = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        assert_eq!(count_table, count_table_old);
    }

    #[test]
    fn test_count_table_by_vs_old_8() {
        let (gsnv1, gsnv2, _) = setup_test_8();
        let count_table = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        let count_table_old = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        assert_eq!(count_table, count_table_old);
    }

    #[test]
    fn test_count_table_by_vs_old_9_nopad() {
        let (gsnv1, gsnv2, _) = setup_test_9_nopad();
        let count_table = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        let count_table_old = count_table_by(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv());
        assert_eq!(count_table, count_table_old);
    }

    #[test]
    fn test_count_table_by_phe_vs_naive_8_nopad() {
        let (gsnv1, gsnv2, phe) = setup_test_8();
        let count_table = count_table_by_phe(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv(), &phe);
        let count_table_naive =
            count_table_by_phe_naive(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv(), &phe);
        assert_eq!(count_table, count_table_naive);
        //assert_eq!(count_table.0, count_table_naive.0);
        //assert_eq!(count_table.1, count_table_naive.1);
        //assert_eq!(count_table.2, count_table_naive.2);
    }

    #[test]
    fn test_count_table_by_phe_vs_naive_9_nopad() {
        let (gsnv1, gsnv2, phe) = setup_test_9_nopad();
        let count_table = count_table_by_phe(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv(), &phe);
        let count_table_naive =
            count_table_by_phe_naive(&gsnv1.as_genot_snv(), &gsnv2.as_genot_snv(), &phe);
        assert_eq!(count_table, count_table_naive);
        //assert_eq!(count_table.0, count_table_naive.0);
        //assert_eq!(count_table.1, count_table_naive.1);

        let sum = count_table.0.sum() + count_table.1.sum();
        assert_eq!(sum, gsnv1.n());
    }

    #[test]
    fn test_sum_table3by3() {
        let count_table = ((2, 0, 1), (2, 1, 0), (1, 1, 1));
        let sum = sum_table3by3(count_table);
        assert_eq!(sum, 9);
    }

    #[test]
    fn test_calculate_r2_self() {
        let (gsnv, _) = setup_test_3();

        let r2 = calc_r2(&gsnv.as_genot_snv(), &gsnv.as_genot_snv());
        assert_float_absolute_eq!(r2, 1.0);
        //assert!(is_eq_f64(r2, 1.0, 1e-7));
    }

    #[test]
    fn test_calculate_r2_2() {
        //let x0: [u8; 4] = [0, 1, 1, 2]; mean: 1, var: 0.5
        //let x1: [u8; 4] = [2, 1, 1, 0]; mean: 1, var: 0.5
        let r2_exp = -1.0;

        let (gsnv0, gsnv1) = setup_test_2();
        let r2 = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());

        assert_float_absolute_eq!(r2, r2_exp);
    }

    #[test]
    fn test_calculate_r2_3() {
        let r2_exp = -0.3464;

        let (gsnv0, gsnv1) = setup_test_3();
        let r2 = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());

        assert_float_absolute_eq!(r2, r2_exp, 1e-4);
    }

    #[test]
    fn test_calculate_r2_4() {
        let (gsnv0, gsnv1) = setup_test_4();
        let r2 = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());

        assert_float_absolute_eq!(r2, 1.0);
    }

    #[test]
    fn test_calculate_r2_5() {
        let (gsnv0, gsnv1) = setup_test_5();
        let r2 = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());
        assert!(r2.is_nan());
    }
    #[test]
    fn test_calculate_r2_6() {
        let (gsnv0, gsnv1) = setup_test_6();
        let r2 = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());
        assert!(r2.is_nan());
    }

    #[test]
    fn test_calculate_r2_7() {
        let (gsnv0, gsnv1) = setup_test_7();
        let r2_test7 = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());

        let (gsnv0, gsnv1) = setup_test_7_missing();
        let r2_test7_missing = calc_r2(&gsnv0.as_genot_snv(), &gsnv1.as_genot_snv());

        assert_float_absolute_eq!(r2_test7, r2_test7_missing)
    }
}
