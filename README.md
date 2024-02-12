# Project Name

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

A brief description of your project.

## Table of Contents

- [Project Name](#project-name)
  - [Description](#description)
  - [Table of Contents](#table-of-contents)
  - [Build](#build)
    - [Versionamento](#versionamento)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Contributing](#contributing)
  - [License](#license)

## Build

Para compilar o projeto em um pacote distribuível, confira se todas as bibliotecas estão instaladas ([`requirements.txt`](/requirements.txt)) e execute:

```shell
hatchling build
```

Isso vai criar na pasta `dist/` o arquivo do pacote. Falta configurar o CI e o PyPI.

### Versionamento

O arquivo [`__about__.py`](/__about__.py) contém a versão atual do projeto. **NÃO O EDITE MANUALMENTE**. Use os comandos a seguir para atualizar a versão
```shell
hatchling version # para ver a versão
hatchling version major # atualiza, por exemplo, 1.1.0 -> 2.0.0 (próxima major)
hatchling version minor # atualiza, por exemplo, 1.1.0 -> 1.2.0 (próxima minor)
hatchling version patch # atualiza, por exemplo, 1.1.0 -> 1.1.1 (próxima patch)
hatchling version 2.1.0 # atualiza para a versão 2.1.0 (proibido regredir a versão)
```

Após fazer isso, faça um commit:
```shell
git add .
git commit -m ":bookmark: BUMP VERSION TO <versao>"
```

## Installation

Instructions on how to install and run your project.

## Usage

Instructions on how to use your project.

## Contributing

Guidelines on how to contribute to your project.

## License

This project is licensed under the [MIT License](LICENSE).
