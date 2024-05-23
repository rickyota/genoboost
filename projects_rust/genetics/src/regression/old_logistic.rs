// [ref](https://github.com/VolodymyrOrlov/vlad-orlov.com_code/blob/main/classifier-in-rust/src/main.rs)

use argmin::prelude::*;
use argmin::solver::linesearch::{ArmijoCondition, BacktrackingLineSearch};
use argmin::solver::quasinewton::LBFGS;
use nalgebra::{DMatrix, DVector};


struct BinaryObjectiveFunction<'a> {
    x: &'a DMatrix<f64>,
    y: &'a DVector<f64>
}

impl<'a> ArgminOp for BinaryObjectiveFunction<'a> {
    type Param = Vec<f64>;
    type Output = f64;
    type Hessian = Vec<Vec<f64>>;
    type Jacobian = ();
    type Float = f64;

    fn apply(&self, w: &Self::Param) -> Result<Self::Output, Error> {
        let mut f = 0f64;
        let (n, _) = self.x.shape();

        for i in 0..n {
            let wx = dot(w, &self.x, i);
            f += self.y[i] * sigmoid(wx).ln() + (1.0 - self.y[i]) * (1.0 - sigmoid(wx)).ln();            
        }
        
        Ok(-f)
    }

    fn gradient(&self, w: &Self::Param) -> Result<Self::Param, Error> {
        let (n, p) = self.x.shape();
        let mut g = vec![0f64; w.len()];

        for i in 0..n {
            let wx = dot(w, &self.x, i);

            let dyi = sigmoid(wx) - self.y[i];
            for j in 0..p {
                g[j] += dyi * self.x[(i, j)];
            }
            g[p] += dyi;
        }
        Ok(g)
    }    
}

fn sigmoid(v: f64) -> f64 {
    if v < -40. {
        0.
    } else if v > 40. {
        1.
    } else {
        1. / (1. + f64::exp(-v))
    }
}

fn dot(w: &Vec<f64>, x: &DMatrix<f64>, m_row: usize) -> f64 {
    let mut sum = 0f64;
    let (_, p) = x.shape();
    for i in 0..p {
        sum += x[(m_row, i)] * w[i];
    }

    sum + w[p]
}

fn optimize(x: &DMatrix<f64>, y: &DVector<f64>) -> Result<(DVector<f64>, f64), Error> {      

    let (_, p) = x.shape();

    // Define cost function
    let cost = BinaryObjectiveFunction { x, y };

    // Define initial parameter vector
    let init_param: Vec<f64> = vec![0f64; p + 1];
    
    // Set condition
    let cond = ArmijoCondition::new(0.5)?;

    // set up a line search
    let linesearch = BacktrackingLineSearch::new(cond).rho(0.9)?;

    // Set up solver
    let solver = LBFGS::new(linesearch, 7);

    // Run solver
    let res = Executor::new(cost, solver, init_param)
        // .add_observer(ArgminSlogLogger::term(), ObserverMode::Always)
        .max_iters(100000)
        .run()?;

    let w = DVector::from_row_slice(&res.state().best_param);
        
    //Ok((w.rows(0, 30).into_owned(), w[30]))
    Ok((w.rows(0, p).into_owned(), w[p]))
}


pub fn logreg(x: Vec<f64>, y:Vec<f64>, row_n:usize, col_n:usize) ->(f64,Vec<f64>){
	println!("create x,y");

	let x=DMatrix::from_vec(row_n,col_n,x);
	let y=DVector::from_vec(y);

    //let file = File::open("creditcard.csv").unwrap();
    //let credit: DMatrix<f64> = parse_csv(BufReader::new(file)).unwrap(); 
    
    //let x = credit.columns(0, 30).into_owned();
    //let y = credit.column(30).into_owned(); 

    //let (x_train, x_test, y_train, y_test) = train_test_split(&x, &y.transpose(), 0.2, true);

	println!("start optimize");

    let (coeff, intercept) = optimize(&x, &y).unwrap();
    //let (coeff, intercept) = optimize(&x, &y.transpose()).unwrap();
    //let (coeff, intercept) = optimize(&x_train, &y_train.transpose()).unwrap();

	println!("done");

    println!("{:?}", coeff);
    println!("{}", intercept);

    //let y_hat = predict(&x_test, &coeff, intercept);
    let coeff = coeff.iter().map(|x| *x).collect::<Vec<f64>>();

    //println!("{}", roc_auc_score(&y_test, &y_hat.transpose()));

	(intercept, coeff)
}