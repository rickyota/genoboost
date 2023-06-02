// supress all warnings in test_rust
#![allow(warnings)]

// https://pyo3.rs/latest/python_from_rust.html

//use pyo3;
use pyo3::prelude::*;
//use pyo3::types::PyList;
//use std::env;

use ndarray;
use ndarray::{array,Array1};


use numpy::{IntoPyArray,PyArray, PyArray1, PyArray2, PyArrayDyn, PyReadonlyArrayDyn};

fn fibonacci_py(n: i64) -> PyResult<i64> {
    Python::with_gil(|py| {
        let fun = PyModule::from_code(
            py,
            r#"
def fibonacci(n):
    a, b = 1, 0
    for _ in range(n):
        a, b = b, a + b
    return b
            "#,
            "",
            "",
        )?
        .getattr("fibonacci")?;

        fun.call1((n,))?.extract()
    })
}

fn fibonacci2_noarg_py() -> PyResult<()> {
    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/py/logreg.py"));
    let from_python = Python::with_gil(|py| -> PyResult<Py<PyAny>> {
        PyModule::from_code(py, py_app, "logreg", "logreg")?;
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "")?
            .getattr("fib5")?
            .into();
        app.call0(py)
    });

    println!("py: {}", from_python?);
    Ok(())
}

// ok!
// [ref](https://docs.rs/pyo3/latest/pyo3/prelude/struct.PyModule.html#method.from_code)
fn fibonacci3_noarg_py() -> i64 {
    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/py/logreg.py"));

    Python::with_gil(|py| -> i64 {
        PyModule::from_code(py, py_app, "logreg", "logreg").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "").unwrap()
            .getattr("fib5").unwrap()
            .into();
        let res = app.call0(py);


        let ans: i64 =res.unwrap().as_ref(py).extract().unwrap();

        //let ans_pyany:&PyAny=res.unwrap().as_ref(py);
        //let ans: i64=ans_pyany.extract().unwrap();
        ans


        // [ref](https://pyo3.rs/latest/types.html)
        // [ref](https://pyo3.rs/latest/types.html#pyt-and-pyobject)
        /*
        let ans_pypyany: Py<PyAny> = res.unwrap();
        let ans_pyany: &PyAny = ans_pypyany.as_ref(py);
        let ans: i64 = ans_pyany.extract().unwrap();
        ans
         */
    })

}


// ok!
// [ref](https://docs.rs/pyo3/latest/pyo3/prelude/struct.PyModule.html#method.from_code)
fn fibonacci3_py(n: i64) -> i64 {
    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/py/logreg.py"));

    Python::with_gil(|py| -> i64 {
        PyModule::from_code(py, py_app, "logreg", "logreg").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "").unwrap()
            .getattr("fib").unwrap()
            .into();
        //let res = app.call0(py);
        let res = app.call1(py,(n,));


        let ans: i64 =res.unwrap().as_ref(py).extract().unwrap();

        ans
    })
}



fn array_sum() -> f64 {

    let v: Array1<f64>=array![1.0,1.2,-0.3];

    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/py/logreg.py"));

    Python::with_gil(|py| -> f64 {
        PyModule::from_code(py, py_app, "logreg", "logreg").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "").unwrap()
            .getattr("ar_sum").unwrap()
            .into();
        //let res = app.call0(py);
        let res = app.call1(py,(v.into_pyarray(py) ,));
        //let res = app.call1(py,(v ,));


        let ans: f64 =res.unwrap().as_ref(py).extract().unwrap();

        ans
    })
}



fn array_twice() -> Vec<f64>{

    //let v: Array1<f64>=array![1.0,1.2,-0.3];
    let v: Vec<f64>=vec![1.0,1.2,-0.3];

    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/py/logreg.py"));

    //Python::with_gil(|py| -> PyArray1<f64> {
    Python::with_gil(|py| -> Vec<f64>{
        PyModule::from_code(py, py_app, "logreg", "logreg").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "").unwrap()
            .getattr("ar_twice").unwrap()
            .into();
        //let res = app.call0(py);
        let res = app.call1(py,(v.into_pyarray(py) ,));
        //let res = app.call1(py,(v ,));


        //let ans=res.unwrap().as_ref(py);
        //let ans: Array1<f64> =res.unwrap().as_ref(py).extract().unwrap();
        let ans: Vec<f64> =res.unwrap().as_ref(py).extract().unwrap();
        //let ans: PyArray1<f64> =res.unwrap().as_ref(py);

        ans
    })
}


fn logreg()->Vec<f64>{

    let vs: Vec<Vec<f64>>=vec![vec![1.0,1.2,-0.3],vec![-1.5,2.1,0.1]];
    let y: Vec<bool>=vec![true,false];
    let w: Vec<f64>=vec![0.9,0.1];

    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/py/logreg.py"));

    //Python::with_gil(|py| -> PyArray1<f64> {
    Python::with_gil(|py| -> Vec<f64>{
        PyModule::from_code(py, py_app, "logreg", "logreg").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "").unwrap()
            .getattr("logreg").unwrap()
            .into();
        //let res = app.call0(py);
        // [ref](https://docs.rs/numpy/latest/numpy/array/struct.PyArray.html#method.from_vec2)
        let args=(PyArray::from_vec2(py, &vs).unwrap(), y.into_pyarray(py),w.into_pyarray(py));
        let res = app.call1(py,args);
        //let res = app.call1(py,(v ,));


        //let ans=res.unwrap().as_ref(py);
        //let ans: Array1<f64> =res.unwrap().as_ref(py).extract().unwrap();
        let ans: Vec<f64> =res.unwrap().as_ref(py).extract().unwrap();
        //let ans: PyArray1<f64> =res.unwrap().as_ref(py);

        ans
    })


}

pub fn test() {
    println!("\ntest");

    //fibonacci2_noarg_py();

    //let x=fibonacci_py(5).unwrap();
    //let x=fibonacci2_py(5);
    //let x=fibonacci2_py(5).unwrap();
    //let x=fibonacci3_noarg_py();
    //let x=fibonacci3_py(7);
    //let x=array_sum();
    //let x=array_twice();
    let x=logreg();
    println!("ans {:?}",x);
}
