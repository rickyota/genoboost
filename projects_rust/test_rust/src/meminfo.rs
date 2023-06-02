use sysinfo::RefreshKind;
use sysinfo::SystemExt;

pub fn test() {
    //println!("create new all");
    // hang up
    //let system = sysinfo::System::new_all();

    println!("create memory");
    // -> solved!!
    // only load memory
    let system = sysinfo::System::new_with_specifics(RefreshKind::new().with_memory());
    println!("created sys");

    let mem = system.available_memory() as usize;
    println!("mem {}", mem);
}
