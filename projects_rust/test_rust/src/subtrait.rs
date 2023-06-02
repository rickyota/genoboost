trait Animal {
    fn get_animal_name(&self) -> String;
}

trait Mammals: Animal {
    //fn get_mammal_name(&self) -> String;
    fn get_mammal_name(&self) -> String {
        self.get_animal_name() + " Default"
    }
}

struct Human {}

impl Animal for Human {
    fn get_animal_name(&self) -> String {
        "Animal".to_owned()
    }
}

impl Mammals for Human {
    fn get_mammal_name(&self) -> String {
        // This is great!!
        self.get_animal_name() + " Mammal"
    }
}

struct Cat {}

impl Animal for Cat {
    fn get_animal_name(&self) -> String {
        "Animal in cat".to_owned()
    }
}

impl Mammals for Cat {
    // use default
    //fn get_mammal_name(&self) -> String {
    //    // This is great!!
    //    self.get_animal_name() + " Mammal"
    //}
}

pub fn test() {
    let a = Human {};
    println!("{}", a.get_mammal_name());

    let b = Cat {};
    println!("{}", b.get_mammal_name());
}
