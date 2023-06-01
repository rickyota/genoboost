// supress all warnings in test_rust
#![allow(warnings)]

mod bufwriter_std;
mod genotype_plan1;
mod genotype_plan2;
mod genotype_plan3;
mod genotype_plan4;
mod genotype_trait;
mod iterator;
mod iterator2;
mod lifetime;
mod logs;
mod meminfo;
mod mut_ref_in_struct;
mod option_asref;
mod path;
mod simd;
mod sort_mut;
mod str_lifetime;
mod struct_asref;
mod subtrait;
mod trait_basic;
mod trait_basic2;
mod trait_boundary;
mod trait_boundary_ok;
mod trait_boundary_ok2;
mod iterator_max_by;
mod trait_boundary_ok_simple;
mod trait_boundary_simple;
mod vec;
mod chain_type;

pub fn test() {
    /*
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
    println!("\ntest_simd");
    simd::test();
    }

    println!("\ntest_predict");
    mut_ref_in_struct::test();

    println!("\ntrait basic");
    trait_basic::test();

    println!("\ntrait basic2");
    trait_basic2::test();

    println!("\ntrait boundary_simple");
    trait_boundary_simple::test();

    println!("\ntrait boundary ok");
    trait_boundary_ok::test();

    println!("\ntrait boundary ok 2");
    trait_boundary_ok2::test();

    println!("\ntrait boundary ok simple");
    trait_boundary_ok_simple::test();

    println!("\ntrait boundary");
    trait_boundary::test();

    println!("\nstruct asref");
    struct_asref::test();

    println!("\noption as ref");
    option_asref::test();

    println!("\nsubtrait");
    subtrait::test();

    println!("\nstr_lifetime");
    str_lifetime::test();

    println!("\nsort mut");
    sort_mut::test();
    */

    /*     println!("\nbufwriter std");
       bufwriter_std::test();
    */
    /*
    println!("\niterator");
    iterator::test();

    println!("\niterator2");
    iterator2::test();

    println!("\ngenotype_plan1");
    genotype_plan1::test();

    println!("\ngenotype_plan2");
    genotype_plan2::test();

    println!("\ngenotype_plan3");
    genotype_plan3::test();

    println!("\ngenotype_plan4");
    genotype_plan4::test();

    println!("\nlifetime");
    lifetime::test();

    println!("\nvec");
    vec::test();

    println!("\ngenotype trait");
    genotype_trait::test();

    println!("\ncreate dir");
    path::test();
    */

    /*     println!("\nlog");
    logs::test(); */

    /*
    println!("\nmeminfo");
    meminfo::test();
     */


    /*
    println!("\nchain_type");
    chain_type::test();
     */


    println!("\niterator_max_by");
    iterator_max_by::test();
    
}
