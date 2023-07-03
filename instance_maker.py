import argparse

from MIP_MCPP.instance import Instance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "name", help="Instance Name, should be '[grid width]x[grid height]-[Characteristics]-k[# of roots]'")
    parser.add_argument("--map", default=None,
                        help="Path to the Binary Map to Create the Instance")

    args = parser.parse_args()

    if args.map:
        istc = Instance.create_from_binary_map(args.map)
    else:
        istc = Instance.create_random_free(args.name)

    if istc:
        istc.write(args.name + '.istc')


if __name__ == "__main__":
    main()
