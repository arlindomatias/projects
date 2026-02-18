CREATE DATABASE laboratory_db
  DEFAULT CHARACTER SET utf8mb4
  COLLATE utf8mb4_unicode_ci;

USE laboratory_db;

CREATE TABLE roles (
  role_id INT PRIMARY KEY,
  role_name VARCHAR(100) NOT NULL
) ENGINE=InnoDB;

CREATE TABLE users (
  user_id INT PRIMARY KEY,
  user_name VARCHAR(150) NOT NULL,
  role_id INT,
  birthday DATE,

  CONSTRAINT fk_users_role
    FOREIGN KEY (role_id)
    REFERENCES roles(role_id)
) ENGINE=InnoDB;

CREATE TABLE room (
  room_id INT PRIMARY KEY,
  laboratory_technician_id INT,
  last_clean TIMESTAMP,

  CONSTRAINT fk_room_technician
    FOREIGN KEY (laboratory_technician_id)
    REFERENCES users(user_id)
) ENGINE=InnoDB;

CREATE TABLE supplier (
  supplier_id INT PRIMARY KEY,
  supplier_name VARCHAR(150) NOT NULL,
  contact_technician_id INT,

  CONSTRAINT fk_supplier_contact
    FOREIGN KEY (contact_technician_id)
    REFERENCES users(user_id)
) ENGINE=InnoDB;

CREATE TABLE experiments (
  experiment_id INT PRIMARY KEY,
  experiment_name VARCHAR(200) NOT NULL,
  room_id INT,

  CONSTRAINT fk_experiment_room
    FOREIGN KEY (room_id)
    REFERENCES room(room_id)
) ENGINE=InnoDB;

CREATE TABLE equipment (
  equipment_id INT PRIMARY KEY,
  equipment_name VARCHAR(200) NOT NULL,
  responsible_technician_id INT,
  manufacturer_id INT,
  last_maintenance TIMESTAMP,

  CONSTRAINT fk_equipment_technician
    FOREIGN KEY (responsible_technician_id)
    REFERENCES users(user_id),

  CONSTRAINT fk_equipment_supplier
    FOREIGN KEY (manufacturer_id)
    REFERENCES supplier(supplier_id)
) ENGINE=InnoDB;

CREATE TABLE experiment_users (
  experiment_id INT,
  user_id INT,

  PRIMARY KEY (experiment_id, user_id),

  CONSTRAINT fk_exp_users_experiment
    FOREIGN KEY (experiment_id)
    REFERENCES experiments(experiment_id)
    ON DELETE CASCADE,

  CONSTRAINT fk_exp_users_user
    FOREIGN KEY (user_id)
    REFERENCES users(user_id)
    ON DELETE CASCADE
) ENGINE=InnoDB;

CREATE TABLE experiment_equipment (
  experiment_id INT,
  equipment_id INT,

  PRIMARY KEY (experiment_id, equipment_id),

  CONSTRAINT fk_exp_equipment_experiment
    FOREIGN KEY (experiment_id)
    REFERENCES experiments(experiment_id)
    ON DELETE CASCADE,

  CONSTRAINT fk_exp_equipment_equipment
    FOREIGN KEY (equipment_id)
    REFERENCES equipment(equipment_id)
    ON DELETE CASCADE
) ENGINE=InnoDB;
